from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Any, Dict, Tuple
import numpy as np

DEFAULT_RTOL = 1e-10
DEFAULT_ATOL = 1e-12
DEFAULT_STABILITY_FACTOR = 10.0
MAX_VERIFIER_RTOL = 1e-8
MAX_VERIFIER_ATOL = 1e-10
NUMERICAL_POLICY_ID = "ATI-numerics-v0.4-rc2"


def validate_tolerances(rtol: float = DEFAULT_RTOL, atol: float = DEFAULT_ATOL,
                        stability_factor: float = DEFAULT_STABILITY_FACTOR,
                        *, verifier: bool = False) -> Tuple[float, float, float]:
    rtol = float(rtol); atol = float(atol); stability_factor = float(stability_factor)
    if not np.isfinite(rtol) or rtol <= 0:
        raise ValueError("rtol must be finite and > 0")
    if not np.isfinite(atol) or atol < 0:
        raise ValueError("atol must be finite and >= 0")
    if not np.isfinite(stability_factor) or stability_factor <= 1:
        raise ValueError("stability_factor must be finite and > 1")
    if verifier:
        if rtol > MAX_VERIFIER_RTOL:
            raise ValueError(f"verifier rtol must be <= {MAX_VERIFIER_RTOL:g}")
        if atol > MAX_VERIFIER_ATOL:
            raise ValueError(f"verifier atol must be <= {MAX_VERIFIER_ATOL:g}")
    return rtol, atol, stability_factor


def absrel_threshold(scale: float, *, rtol: float = DEFAULT_RTOL, atol: float = DEFAULT_ATOL) -> float:
    rtol, atol, _ = validate_tolerances(rtol, atol)
    scale = float(scale)
    if not np.isfinite(scale) or scale < 0:
        raise ValueError("scale must be finite and nonnegative")
    return float(atol + rtol * max(scale, 1.0))


def rank_threshold(singular_values, shape, *, rtol: float = DEFAULT_RTOL, atol: float = DEFAULT_ATOL) -> float:
    rtol, atol, _ = validate_tolerances(rtol, atol)
    s = np.asarray(singular_values, dtype=np.float64)
    smax = float(s[0]) if s.size else 0.0
    return float(atol + rtol * max(shape) * max(smax, 1.0))


def tight_loose(base_threshold: float, stability_factor: float = DEFAULT_STABILITY_FACTOR):
    _, _, stability_factor = validate_tolerances(DEFAULT_RTOL, DEFAULT_ATOL, stability_factor)
    t = float(base_threshold)
    if not np.isfinite(t) or t < 0:
        raise ValueError("base_threshold must be finite and nonnegative")
    return t / stability_factor, t * stability_factor


def classify_nonnegative_residual(value: float, base_threshold: float,
                                  stability_factor: float = DEFAULT_STABILITY_FACTOR):
    value = float(value)
    if value < 0 or not np.isfinite(value):
        raise ValueError("residual must be finite and nonnegative")
    tight, loose = tight_loose(base_threshold, stability_factor)
    if value <= tight:
        return "ZERO_WITHIN_NUMERICAL_MODEL", tight, loose
    if value > loose:
        return "NONZERO_WITHIN_NUMERICAL_MODEL", tight, loose
    return "UNKNOWN_NUMERICAL", tight, loose


def svd_rank_info(A, *, rtol: float = DEFAULT_RTOL, atol: float = DEFAULT_ATOL,
                  stability_factor: float = DEFAULT_STABILITY_FACTOR,
                  full_matrices: bool = False) -> Dict[str, Any]:
    rtol, atol, stability_factor = validate_tolerances(rtol, atol, stability_factor)
    A = np.asarray(A, dtype=np.float64)
    if A.ndim != 2 or not np.all(np.isfinite(A)):
        raise ValueError("A must be a finite matrix")
    u, s, vh = np.linalg.svd(A, full_matrices=full_matrices)
    base = rank_threshold(s, A.shape, rtol=rtol, atol=atol)
    tight, loose = tight_loose(base, stability_factor)
    rank_tight = int(np.sum(s > tight))
    rank_loose = int(np.sum(s > loose))
    stable = rank_tight == rank_loose
    rank = rank_loose if stable else None
    return {
        "u": u, "s": s, "vh": vh,
        "rank": rank,
        "rank_tight": rank_tight,
        "rank_loose": rank_loose,
        "rank_stable": bool(stable),
        "rank_threshold_base": float(base),
        "rank_threshold_tight": float(tight),
        "rank_threshold_loose": float(loose),
    }


def pinv_from_svd(A, *, rtol: float = DEFAULT_RTOL, atol: float = DEFAULT_ATOL,
                  stability_factor: float = DEFAULT_STABILITY_FACTOR,
                  require_stable_rank: bool = True):
    info = svd_rank_info(A, rtol=rtol, atol=atol, stability_factor=stability_factor, full_matrices=False)
    if require_stable_rank and not info["rank_stable"]:
        return None, info
    rank = info["rank_loose"] if info["rank"] is None else info["rank"]
    u, s, vh = info["u"], info["s"], info["vh"]
    invs = np.zeros_like(s)
    if rank:
        invs[:rank] = 1.0 / s[:rank]
    pinv = (vh.T * invs) @ u.T
    return pinv, info


def policy_record(*, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                  stability_factor=DEFAULT_STABILITY_FACTOR):
    rtol, atol, stability_factor = validate_tolerances(rtol, atol, stability_factor)
    return {
        "policy_id": NUMERICAL_POLICY_ID,
        "rtol": rtol,
        "atol": atol,
        "stability_factor": stability_factor,
        "rank_threshold": "atol + rtol*max(shape)*max(smax,1)",
        "zero_nonzero_decision": "tight=base/stability_factor; loose=base*stability_factor; middle=UNKNOWN",
    }
