from __future__ import annotations

import hashlib
from typing import Iterable, Dict, Any, Optional, Mapping
import numpy as np

from .numerics import (
    DEFAULT_RTOL, DEFAULT_ATOL, DEFAULT_STABILITY_FACTOR,
    validate_tolerances, absrel_threshold, rank_threshold, tight_loose,
    svd_rank_info, pinv_from_svd, policy_record,
)

_EPS = np.finfo(float).eps


def _as_matrix(X, name="operator"):
    X = np.asarray(X, dtype=np.float64)
    if X.ndim != 2:
        raise ValueError(f"{name} must be a 2D matrix")
    if not np.all(np.isfinite(X)):
        raise ValueError(f"{name} contains non-finite values")
    return np.ascontiguousarray(X)


def validate_operator(X, name="operator") -> Dict[str, Any]:
    X = _as_matrix(X, name)
    row_norms = np.linalg.norm(X, axis=1)
    return {
        "name": name, "shape": list(X.shape), "finite": True,
        "zero_rows": np.flatnonzero(row_norms == 0).astype(int).tolist(),
        "row_norm_min": float(row_norms.min()) if row_norms.size else 0.0,
        "row_norm_max": float(row_norms.max()) if row_norms.size else 0.0,
        "operator_sha256": array_sha256(X),
    }


def array_sha256(X) -> str:
    X = _as_matrix(X)
    h = hashlib.sha256()
    h.update(str(X.shape).encode("ascii")); h.update(b"|float64|C|")
    h.update(X.astype("<f8", copy=False).tobytes(order="C"))
    return h.hexdigest()


def _rank_threshold(s, shape, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    return rank_threshold(s, shape, rtol=rtol, atol=atol)


def numerical_rank(A, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                   stability_factor=DEFAULT_STABILITY_FACTOR):
    A = _as_matrix(A, "A")
    info = svd_rank_info(A, rtol=rtol, atol=atol, stability_factor=stability_factor)
    return int(info["rank_loose"] if info["rank"] is None else info["rank"])


def rowspace_basis(A, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                   stability_factor=DEFAULT_STABILITY_FACTOR):
    A = _as_matrix(A, "A")
    info = svd_rank_info(A, rtol=rtol, atol=atol, stability_factor=stability_factor, full_matrices=False)
    r = info["rank_loose"] if info["rank"] is None else info["rank"]
    return info["vh"][:r, :], info["s"], info["rank_threshold_base"]


def rowspace_projector(A, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                       stability_factor=DEFAULT_STABILITY_FACTOR):
    A = _as_matrix(A, "A")
    info = svd_rank_info(A, rtol=rtol, atol=atol, stability_factor=stability_factor, full_matrices=False)
    if not info["rank_stable"]:
        return None, info
    r = info["rank"]
    Q = info["vh"][:r, :]
    n = A.shape[1]
    P = np.zeros((n, n), dtype=np.float64) if r == 0 else Q.T @ Q
    return P, info


def principal_angles(A, B, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                     stability_factor=DEFAULT_STABILITY_FACTOR):
    A = _as_matrix(A, "A"); B = _as_matrix(B, "B")
    if A.shape[1] != B.shape[1]:
        raise ValueError("A and B must act on the same latent dimension")
    ia = svd_rank_info(A, rtol=rtol, atol=atol, stability_factor=stability_factor, full_matrices=False)
    ib = svd_rank_info(B, rtol=rtol, atol=atol, stability_factor=stability_factor, full_matrices=False)
    if not ia["rank_stable"] or not ib["rank_stable"]:
        return np.array([], dtype=np.float64)
    QA = ia["vh"][:ia["rank"], :]; QB = ib["vh"][:ib["rank"], :]
    if QA.shape[0] == 0 or QB.shape[0] == 0:
        return np.array([], dtype=np.float64)
    c = np.linalg.svd(QA @ QB.T, compute_uv=False)
    return np.arccos(np.clip(c, 0.0, 1.0))


def _relative_residual(B, P):
    if P is None:
        return None, None, None, None
    R = B @ (np.eye(B.shape[1]) - P)
    nr = float(np.linalg.norm(R, ord=2)) if R.size else 0.0
    nb = float(np.linalg.norm(B, ord=2)) if B.size else 0.0
    rel = nr / max(nb, _EPS)
    return R, nr, nb, rel


def interoperability_deficit(A, B, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                             stability_factor=DEFAULT_STABILITY_FACTOR) -> Dict[str, Any]:
    A = _as_matrix(A, "A"); B = _as_matrix(B, "B")
    if A.shape[1] != B.shape[1]: raise ValueError("A and B must act on the same latent dimension")
    P, ia = rowspace_projector(A, rtol, atol, stability_factor)
    if P is None:
        return {"status": "UNKNOWN_NUMERICAL", "reason": "source_rank_unstable",
                "missing_dimension_numerical": None, "residual_relative": None,
                "principal_angles_radians": np.array([], dtype=np.float64)}
    R, nr, nb, rel = _relative_residual(B, P)
    ir = svd_rank_info(R, rtol=rtol, atol=atol, stability_factor=stability_factor, full_matrices=False)
    return {
        "status": "STABLE" if ir["rank_stable"] else "UNKNOWN_NUMERICAL",
        "residual_norm2": nr, "residual_relative": rel,
        "missing_dimension_numerical": int(ir["rank"]) if ir["rank_stable"] else None,
        "residual_singular_values": ir["s"],
        "principal_angles_radians": principal_angles(A, B, rtol, atol, stability_factor),
        "residual_rank_stable": ir["rank_stable"],
    }


def find_kernel_witness(A, B, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                        stability_factor=DEFAULT_STABILITY_FACTOR) -> Optional[Dict[str, Any]]:
    A = _as_matrix(A, "A"); B = _as_matrix(B, "B")
    if A.shape[1] != B.shape[1]: raise ValueError("A and B must act on the same latent dimension")
    info = svd_rank_info(A, rtol=rtol, atol=atol, stability_factor=stability_factor, full_matrices=True)
    if not info["rank_stable"]: return None
    r = info["rank"]
    N = info["vh"][r:, :].T
    if N.shape[1] == 0: return None
    BN = B @ N
    if BN.size == 0: return None
    _, sb, vhb = np.linalg.svd(BN, full_matrices=False)
    if sb.size == 0: return None
    bscale = max(float(np.linalg.norm(B, 2)), 1.0)
    base = absrel_threshold(bscale, rtol=rtol, atol=atol)
    tight, loose = tight_loose(base, stability_factor)
    if float(sb[0]) <= loose: return None
    v = N @ vhb[0, :]
    nv = float(np.linalg.norm(v))
    if not np.isfinite(nv) or nv == 0: return None
    v = v / nv
    return {
        "v": v, "A_v_norm2": float(np.linalg.norm(A @ v)),
        "B_v_norm2": float(np.linalg.norm(B @ v)), "gain": float(sb[0]),
        "nullity_A": int(N.shape[1]), "rank_A": int(r),
        "rank_threshold_A": float(info["rank_threshold_base"]),
    }


def exact_recovery_analysis(A, B, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                            stability_factor=DEFAULT_STABILITY_FACTOR) -> Dict[str, Any]:
    """Fail-closed numerical row-space containment analysis."""
    rtol, atol, stability_factor = validate_tolerances(rtol, atol, stability_factor)
    A = _as_matrix(A, "A"); B = _as_matrix(B, "B")
    if A.shape[1] != B.shape[1]: raise ValueError("A and B must act on the same latent dimension")

    P, ia = rowspace_projector(A, rtol, atol, stability_factor)
    ib = svd_rank_info(B, rtol=rtol, atol=atol, stability_factor=stability_factor, full_matrices=False)
    nb = float(np.linalg.norm(B, 2)) if B.size else 0.0
    baseB = absrel_threshold(nb, rtol=rtol, atol=atol)
    tightB, looseB = tight_loose(baseB, stability_factor)

    if P is None or not ib["rank_stable"]:
        return {
            "numerical_status":"UNKNOWN_NUMERICAL", "ati_class_unrestricted_psd":"UNKNOWN_NUMERICAL",
            "rank_A": ia["rank_loose"], "rank_B": ib["rank_loose"],
            "rank_A_stable": bool(ia["rank_stable"]), "rank_B_stable": bool(ib["rank_stable"]),
            "rowspace_residual_norm2": None, "rowspace_residual_relative": None,
            "factorization_residual_norm2": None, "factorization_residual_relative": None,
            "M": None, "R": None, "witness": None,
            "numerical_policy": policy_record(rtol=rtol, atol=atol, stability_factor=stability_factor),
            "warning":"Numerical rank is unstable; no exact/no-go claim is permitted."
        }

    R, nr, _, rel = _relative_residual(B, P)
    if nr <= tightB:
        status = "EXACT_WITHIN_NUMERICAL_MODEL"
    elif nr > looseB:
        status = "NONEXACT_WITHIN_NUMERICAL_MODEL"
    else:
        status = "UNKNOWN_NUMERICAL"

    pinvA, _ = pinv_from_svd(A, rtol=rtol, atol=atol, stability_factor=stability_factor, require_stable_rank=True)
    M = None if pinvA is None else B @ pinvA
    factor_norm = factor_rel = None
    if M is not None:
        FR = B - M @ A
        factor_norm = float(np.linalg.norm(FR, 2)) if FR.size else 0.0
        factor_rel = factor_norm / max(nb, _EPS)
        if status == "EXACT_WITHIN_NUMERICAL_MODEL" and factor_norm > looseB:
            status = "UNKNOWN_NUMERICAL"

    witness = find_kernel_witness(A, B, rtol, atol, stability_factor) if status != "EXACT_WITHIN_NUMERICAL_MODEL" else None
    if status == "NONEXACT_WITHIN_NUMERICAL_MODEL" and witness is None:
        status = "UNKNOWN_NUMERICAL"

    ati = ("ATI0_EXACT_ALGEBRAIC_NUMERICALLY_SUPPORTED" if status == "EXACT_WITHIN_NUMERICAL_MODEL"
           else "ATI5_UNBOUNDED_UNDER_UNRESTRICTED_PSD_K" if status == "NONEXACT_WITHIN_NUMERICAL_MODEL"
           else "UNKNOWN_NUMERICAL")
    return {
        "numerical_status": status, "ati_class_unrestricted_psd": ati,
        "rank_A": int(ia["rank"]), "rank_B": int(ib["rank"]),
        "rank_A_stable": True, "rank_B_stable": True,
        "rowspace_residual_norm2": nr, "rowspace_residual_relative": rel,
        "factorization_residual_norm2": factor_norm, "factorization_residual_relative": factor_rel,
        "M": M, "R": R, "witness": witness,
        "residual_threshold_tight": tightB, "residual_threshold_loose": looseB,
        "numerical_policy": policy_record(rtol=rtol, atol=atol, stability_factor=stability_factor),
        "principal_angles_radians": principal_angles(A, B, rtol, atol, stability_factor),
        "warning":"Floating-point classification. Formal symbolic exactness requires exact/certified arithmetic."
    }


def classify_transfer(A, B, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                      stability_factor=DEFAULT_STABILITY_FACTOR):
    return exact_recovery_analysis(A, B, rtol, atol, stability_factor)


def exact_recovery_map(A, B, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    an = exact_recovery_analysis(A, B, rtol, atol)
    if an["numerical_status"] != "EXACT_WITHIN_NUMERICAL_MODEL":
        raise ValueError(f"transfer is not certified exact: {an['numerical_status']}")
    return an["M"]


def trace_bounded_projection_error_bound(A, B, tau, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    if not np.isfinite(tau) or tau < 0: raise ValueError("tau must be finite and nonnegative")
    A = _as_matrix(A,"A"); B=_as_matrix(B,"B")
    if A.shape[1] != B.shape[1]: raise ValueError("A and B must act on same latent dimension")
    P, info = rowspace_projector(A, rtol, atol)
    if P is None: return {"status":"UNKNOWN_NUMERICAL", "bound":None, "reason":"source_rank_unstable"}
    BP=B@P; R=B@(np.eye(B.shape[1])-P)
    bp=float(np.linalg.norm(BP,2)) if BP.size else 0.0; rr=float(np.linalg.norm(R,2)) if R.size else 0.0
    return {"status":"BOUNDED_CONDITIONAL", "tau":float(tau), "BP_norm2":bp, "R_norm2":rr,
            "bound":float(tau)*(2*bp*rr+rr*rr)}


def _packed_sym_values(n: int): return int(n*(n+1)//2)


def minimal_augmentation(A, B, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL) -> Dict[str, Any]:
    A=_as_matrix(A,"A"); B=_as_matrix(B,"B")
    if A.shape[1]!=B.shape[1]: raise ValueError("A and B must act on same latent dimension")
    P, ia=rowspace_projector(A,rtol,atol)
    if P is None:
        return {"status":"UNKNOWN_NUMERICAL","U_add":None,"minimal_dimension_numerical":None,"reason":"source_rank_unstable"}
    R=B@(np.eye(A.shape[1])-P)
    ir=svd_rank_info(R,rtol=rtol,atol=atol,stability_factor=DEFAULT_STABILITY_FACTOR,full_matrices=False)
    if not ir["rank_stable"]:
        return {"status":"UNKNOWN_NUMERICAL","U_add":None,"minimal_dimension_numerical":None,"reason":"residual_rank_unstable"}
    k=int(ir["rank"]); U_add=ir["vh"][:k,:]
    combined=np.vstack([A,U_add]) if k else A.copy()
    check=exact_recovery_analysis(combined,B,rtol,atol)
    a=A.shape[0]
    archive={
        "required_statistic":"full_joint_covariance_of_stacked_[A;U_add]_outputs",
        "separate_C_AA_and_C_UU_sufficient":False,
        "cross_covariance_required":True,
        "source_output_dimension":int(a), "augmentation_dimension":k,
        "incremental_unique_float64_values":int(a*k + _packed_sym_values(k)),
        "incremental_float64_bytes_per_subject":int(8*(a*k + _packed_sym_values(k))),
        "total_joint_unique_float64_values":_packed_sym_values(a+k),
        "total_joint_float64_bytes_per_subject":int(8*_packed_sym_values(a+k)),
        "operator_basis_storage_is_global_not_per_subject":True,
    }
    return {"status":"PASS" if check["numerical_status"]=="EXACT_WITHIN_NUMERICAL_MODEL" else "UNKNOWN_NUMERICAL",
            "U_add":U_add,"minimal_dimension_numerical":k,"rank_residual":k,
            "residual_singular_values":ir["s"],"rank_threshold_residual":ir["rank_threshold_base"],
            "verification_status":check["numerical_status"],"verification_relative_residual":check["rowspace_residual_relative"],
            "archive_requirement":archive}


def design_minimal_augmentation(A,B,rtol=DEFAULT_RTOL,atol=DEFAULT_ATOL): return minimal_augmentation(A,B,rtol,atol)


def universal_representation(operators: Iterable[np.ndarray], rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    ops=[_as_matrix(x,f"operator_{i}") for i,x in enumerate(operators)]
    if not ops: raise ValueError("at least one operator is required")
    n=ops[0].shape[1]
    if any(x.shape[1]!=n for x in ops): raise ValueError("all operators must act on same latent dimension")
    stack=np.vstack(ops)
    info=svd_rank_info(stack,rtol=rtol,atol=atol,stability_factor=DEFAULT_STABILITY_FACTOR,full_matrices=False)
    if not info["rank_stable"]:
        return {"status":"UNKNOWN_NUMERICAL","U":None,"dimension_numerical":None,"reason":"stack_rank_unstable"}
    U=info["vh"][:info["rank"],:]
    checks=[exact_recovery_analysis(U,x,rtol,atol)["numerical_status"] for x in ops]
    d=int(U.shape[0])
    return {"status":"PASS" if all(x=="EXACT_WITHIN_NUMERICAL_MODEL" for x in checks) else "UNKNOWN_NUMERICAL",
            "U":U,"dimension_numerical":d,"stack_rank_numerical":d,"stack_singular_values":info["s"],
            "rank_threshold_stack":info["rank_threshold_base"],"target_verification_statuses":checks,
            "archive_requirement":{
                "for_future_linear_representations_within_joint_span":"store_full_covariance_of_universal_representation",
                "universal_dimension":d,"packed_unique_float64_values":_packed_sym_values(d),
                "packed_float64_bytes_per_subject":int(8*_packed_sym_values(d)),
                "fixed_target_only_note":"If only predeclared target covariance blocks are needed, storing each block can be cheaper than a universal joint representation; universal archive cost is justified by future-span flexibility, not by fixed-target compression."
            }}


def design_universal_representation(operators,rtol=DEFAULT_RTOL,atol=DEFAULT_ATOL): return universal_representation(operators,rtol,atol)


def build_interoperability_graph(operators: Mapping[str,np.ndarray], rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL):
    ids=list(operators); ops={k:_as_matrix(operators[k],k) for k in ids}
    if len({x.shape[1] for x in ops.values()})!=1: raise ValueError("all operators must act on same latent dimension")
    edges=[]; pair_status={}
    for a in ids:
        for b in ids:
            if a==b: continue
            st=exact_recovery_analysis(ops[a],ops[b],rtol,atol)["numerical_status"]
            pair_status[f"{a}->{b}"]=st
            if st=="EXACT_WITHIN_NUMERICAL_MODEL": edges.append([a,b])
    equivalence=[]
    for i,a in enumerate(ids):
        for b in ids[i+1:]:
            if [a,b] in edges and [b,a] in edges: equivalence.append([a,b])
    return {"nodes":ids,"edges":edges,"equivalent_pairs":equivalence,"pair_status":pair_status}


def gram_operator_rank_info(G, operator_shape, *, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                             stability_factor=DEFAULT_STABILITY_FACTOR):
    """Recover operator singular-value/rank information from G = W W^T using W's original shape.

    The original operator shape is mandatory because the shared rank threshold depends on max(shape).
    Gram dimensions alone cannot recover the latent column dimension and therefore cannot reproduce
    the direct numerical contract at large n.
    """
    G=_as_matrix(G,"G")
    if G.shape[0]!=G.shape[1]:
        raise ValueError("G must be square")
    if operator_shape is None or len(tuple(operator_shape))!=2:
        return {"status":"UNKNOWN_NUMERICAL","reason":"operator_shape_required"}
    m,n=(int(operator_shape[0]),int(operator_shape[1]))
    if m!=G.shape[0] or m<0 or n<0:
        raise ValueError("operator_shape is inconsistent with Gram dimensions")
    asym=float(np.linalg.norm(G-G.T,2)) if G.size else 0.0
    scale=float(np.linalg.norm(G,2)) if G.size else 0.0
    asym_tol=absrel_threshold(scale,rtol=rtol,atol=atol)*stability_factor
    if asym>asym_tol:
        return {"status":"UNKNOWN_NUMERICAL","reason":"gram_asymmetry","asymmetry_norm2":asym,"asymmetry_tolerance":asym_tol}
    Gs=0.5*(G+G.T)
    e,Q=np.linalg.eigh(Gs)
    escale=max(float(np.max(np.abs(e))) if e.size else 0.0,1.0)
    negtol=absrel_threshold(escale,rtol=rtol,atol=atol)*stability_factor
    if e.size and float(e.min()) < -negtol:
        return {"status":"UNKNOWN_NUMERICAL","reason":"gram_significant_negative_eigenvalue","min_eigenvalue":float(e.min()),"negative_eigenvalue_tolerance":negtol}
    ep=np.clip(e,0.0,None)
    svals=np.sqrt(ep)[::-1]
    base=rank_threshold(svals,(m,n),rtol=rtol,atol=atol)
    tight,loose=tight_loose(base,stability_factor)
    rank_tight=int(np.sum(svals>tight)); rank_loose=int(np.sum(svals>loose))
    stable=rank_tight==rank_loose
    rank=rank_loose if stable else None
    return {"status":"STABLE" if stable else "UNKNOWN_NUMERICAL","rank":rank,
            "rank_tight":rank_tight,"rank_loose":rank_loose,"rank_stable":bool(stable),
            "singular_values":svals,"rank_threshold_base":float(base),
            "rank_threshold_tight":float(tight),"rank_threshold_loose":float(loose),
            "eigenvalues":e,"eigenvectors":Q,"symmetrized_gram":Gs,
            "operator_shape":[m,n]}


def _gram_pinv_for_operator(G, operator_shape, *, rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                            stability_factor=DEFAULT_STABILITY_FACTOR):
    info=gram_operator_rank_info(G,operator_shape,rtol=rtol,atol=atol,stability_factor=stability_factor)
    if info.get("status")!="STABLE":
        return None,info
    e=np.clip(info["eigenvalues"],0.0,None); Q=info["eigenvectors"]
    # np.linalg.eigh is ascending. Keep exactly the singular directions admitted by the shared
    # direct-operator rank threshold, i.e. sqrt(lambda) > threshold.
    base=info["rank_threshold_base"]
    keep=np.sqrt(e)>base
    inv=np.zeros_like(e); inv[keep]=1.0/e[keep]
    return (Q*inv)@Q.T,info


def gram_containment_analysis(Gaa, Gbb, Gba, *, source_shape=None, target_shape=None,
                              rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                              stability_factor=DEFAULT_STABILITY_FACTOR):
    """Fail-closed containment from Gram blocks under the same numerical contract as direct operators.

    ``source_shape`` and ``target_shape`` are mandatory. A Gram matrix does not encode the latent
    column dimension, while the ATI numerical rank threshold explicitly depends on the original
    operator shape. Omitting shapes therefore returns UNKNOWN rather than silently applying a
    different numerical model.
    """
    rtol,atol,stability_factor=validate_tolerances(rtol,atol,stability_factor)
    Gaa=_as_matrix(Gaa,"Gaa"); Gbb=_as_matrix(Gbb,"Gbb"); Gba=_as_matrix(Gba,"Gba")
    if Gaa.shape[0]!=Gaa.shape[1] or Gbb.shape[0]!=Gbb.shape[1] or Gba.shape!=(Gbb.shape[0],Gaa.shape[0]):
        raise ValueError("invalid Gram block shapes")
    if source_shape is None or target_shape is None:
        return {"numerical_status":"UNKNOWN_NUMERICAL","reason":"original_operator_shapes_required_for_shared_rank_policy"}
    source_shape=tuple(int(x) for x in source_shape); target_shape=tuple(int(x) for x in target_shape)
    if len(source_shape)!=2 or len(target_shape)!=2 or source_shape[0]!=Gaa.shape[0] or target_shape[0]!=Gbb.shape[0] or source_shape[1]!=target_shape[1]:
        raise ValueError("source_shape/target_shape inconsistent with Gram blocks")

    # Exact self-identity control.  When the three Gram blocks are exactly the same
    # and the declared operator shapes agree, the source and target operators are
    # identical in row-space geometry: ||B-A||_F^2 = tr(Gbb)+tr(Gaa)-2tr(Gba)=0.
    # Do not manufacture a tiny nonzero Schur-complement residual by numerically
    # forming G - G G^+ G.  This branch is deliberately exact/structural and does
    # not relax any tolerance for non-identical inputs.
    if source_shape==target_shape and np.array_equal(Gaa,Gbb) and np.array_equal(Gaa,Gba):
        info=gram_operator_rank_info(Gaa,source_shape,rtol=rtol,atol=atol,stability_factor=stability_factor)
        rank=info.get("rank") if info.get("status")=="STABLE" else None
        return {"numerical_status":"EXACT_WITHIN_NUMERICAL_MODEL",
                "reason":"identical_gram_self_identity",
                "direct_residual_norm2":0.0,"target_norm2":float(np.sqrt(max(float(np.linalg.eigvalsh(0.5*(Gbb+Gbb.T)).max()) if Gbb.size else 0.0,0.0))),
                "direct_residual_relative":0.0,"residual_energy_ratio":0.0,
                "residual_gram_min_eigenvalue":0.0,
                "residual_threshold_tight":0.0,"residual_threshold_loose":0.0,
                "source_rank_numerical":rank,"target_rank_numerical":rank,
                "source_rank_stable":bool(info.get("rank_stable",False)),"target_rank_stable":bool(info.get("rank_stable",False)),
                "missing_dimension_numerical":0,"residual_rank_stable":True,
                "source_rank_threshold_base":info.get("rank_threshold_base"),
                "target_rank_threshold_base":info.get("rank_threshold_base"),
                "residual_rank_threshold_base":0.0,
                "source_shape":list(source_shape),"target_shape":list(target_shape),
                "warning":"Exact identical-Gram self identity certified structurally; no Schur-complement subtraction performed."}

    pinv, ia=_gram_pinv_for_operator(Gaa,source_shape,rtol=rtol,atol=atol,stability_factor=stability_factor)
    ib=gram_operator_rank_info(Gbb,target_shape,rtol=rtol,atol=atol,stability_factor=stability_factor)
    if pinv is None or ia.get("status")!="STABLE":
        return {"numerical_status":"UNKNOWN_NUMERICAL","reason":"Gaa_rank_unstable","source_rank_info":ia}
    if ib.get("status")!="STABLE":
        return {"numerical_status":"UNKNOWN_NUMERICAL","reason":"Gbb_rank_unstable","source_rank_info":ia,"target_rank_info":ib}

    raw=Gbb-Gba@pinv@Gba.T
    R=0.5*(raw+raw.T)
    ir=gram_operator_rank_info(R,target_shape,rtol=rtol,atol=atol,stability_factor=stability_factor)
    if ir.get("reason") in {"gram_significant_negative_eigenvalue","gram_asymmetry"}:
        return {"numerical_status":"UNKNOWN_NUMERICAL","reason":"residual_"+ir["reason"],"source_rank_info":ia,"target_rank_info":ib,"residual_rank_info":ir}

    er=np.linalg.eigvalsh(R); eb=np.linalg.eigvalsh(0.5*(Gbb+Gbb.T))
    lamr=max(float(er.max()) if er.size else 0.0,0.0); lamb=max(float(eb.max()) if eb.size else 0.0,0.0)
    nr=float(np.sqrt(lamr)); nb=float(np.sqrt(lamb))
    direct_rel=nr/max(nb,_EPS)
    baseB=absrel_threshold(nb,rtol=rtol,atol=atol); tightB,looseB=tight_loose(baseB,stability_factor)
    status="EXACT_WITHIN_NUMERICAL_MODEL" if nr<=tightB else "NONEXACT_WITHIN_NUMERICAL_MODEL" if nr>looseB else "UNKNOWN_NUMERICAL"
    if not ir.get("rank_stable",False) and status!="EXACT_WITHIN_NUMERICAL_MODEL":
        # Nonzero residual can still be certified from its norm even if its dimension is unstable;
        # preserve the no-go status but abstain on missing-dimension claims.
        missing=None
    else:
        missing=int(ir["rank"]) if ir.get("rank") is not None else None
    return {"numerical_status":status,"direct_residual_norm2":nr,"target_norm2":nb,
            "direct_residual_relative":direct_rel,"residual_energy_ratio":float(lamr/max(lamb,_EPS)),
            "residual_gram_min_eigenvalue":float(er.min()) if er.size else 0.0,
            "residual_threshold_tight":tightB,"residual_threshold_loose":looseB,
            "source_rank_numerical":ia.get("rank"),"target_rank_numerical":ib.get("rank"),
            "source_rank_stable":ia.get("rank_stable"),"target_rank_stable":ib.get("rank_stable"),
            "missing_dimension_numerical":missing,"residual_rank_stable":ir.get("rank_stable",False),
            "source_rank_threshold_base":ia.get("rank_threshold_base"),
            "target_rank_threshold_base":ib.get("rank_threshold_base"),
            "residual_rank_threshold_base":ir.get("rank_threshold_base"),
            "source_shape":list(source_shape),"target_shape":list(target_shape),
            "warning":"Gram classification uses original operator shapes so its rank and residual thresholds match the direct-operator numerical contract. residual_energy_ratio is squared relative spectral residual."}


def tolerance_sweep(A,B,rtols=None,atols=None,stability_factor=DEFAULT_STABILITY_FACTOR):
    if rtols is None: rtols=[1e-12,1e-11,1e-10,1e-9,1e-8]
    if atols is None: atols=[1e-14,1e-13,1e-12,1e-11,1e-10]
    rows=[]; statuses=set()
    for r in rtols:
        for a in atols:
            x=exact_recovery_analysis(A,B,float(r),float(a),stability_factor)
            rows.append({"rtol":float(r),"atol":float(a),"status":x["numerical_status"],"rowspace_residual_relative":x["rowspace_residual_relative"],"rank_A":x["rank_A"],"rank_B":x["rank_B"]})
            statuses.add(x["numerical_status"])
    summary="STABLE_EXACT" if statuses=={"EXACT_WITHIN_NUMERICAL_MODEL"} else "STABLE_NONEXACT" if statuses=={"NONEXACT_WITHIN_NUMERICAL_MODEL"} else "STABLE_UNKNOWN" if statuses=={"UNKNOWN_NUMERICAL"} else "UNSTABLE_ACROSS_TOLERANCE_GRID"
    return {"summary":summary,"statuses":sorted(statuses),"grid":rows}


def generate_certificate(A,B,source_id="source",target_id="target",rtol=DEFAULT_RTOL,atol=DEFAULT_ATOL):
    A=_as_matrix(A,"A"); B=_as_matrix(B,"B"); an=exact_recovery_analysis(A,B,rtol,atol); deficit=interoperability_deficit(A,B,rtol,atol)
    cert={"schema_version":"ATI-certificate-v0.3","producer_version":"ati-reference-0.4.0-rc6",
          "source":{"id":source_id,"shape":list(A.shape),"operator_sha256":array_sha256(A)},
          "target":{"id":target_id,"shape":list(B.shape),"operator_sha256":array_sha256(B)},
          "observation_model":"covariance_congruence_C_equals_W_K_WT","admissible_latent_family":"all_real_symmetric_PSD_K",
          "numerical_policy":an["numerical_policy"],
          "numerics":{"status":an["numerical_status"],"rowspace_residual_norm2":an["rowspace_residual_norm2"],"rowspace_residual_relative":an["rowspace_residual_relative"],"factorization_residual_norm2":an["factorization_residual_norm2"],"factorization_residual_relative":an["factorization_residual_relative"],"rank_A":an["rank_A"],"rank_B":an["rank_B"],"missing_dimension_numerical":deficit.get("missing_dimension_numerical"),"principal_angles_radians":an.get("principal_angles_radians",np.array([])).tolist()},
          "claim":{"ati_class":an["ati_class_unrestricted_psd"],"permitted":[],"prohibited":[]},
          "verification":{"producer_self_check":True,"independent_verifier_required":True,"verifier_must_ignore_producer_tolerances_for_acceptance":True}}
    if an["numerical_status"]=="EXACT_WITHIN_NUMERICAL_MODEL":
        cert["exact_map"]={"M":an["M"].tolist()}; cert["claim"]["permitted"]=["Exact covariance recovery for the declared digital operators within the frozen numerical policy."]; cert["claim"]["prohibited"]=["Formal symbolic exactness unless separately proven.","Correlation-space recovery unless separately established."]
    elif an["numerical_status"]=="NONEXACT_WITHIN_NUMERICAL_MODEL" and an["witness"] is not None:
        cert["no_go_witness"]={"v":an["witness"]["v"].tolist(),"A_v_norm2":an["witness"]["A_v_norm2"],"B_v_norm2":an["witness"]["B_v_norm2"]}; cert["claim"]["permitted"]=["No universal covariance recovery over unrestricted PSD K for the declared digital operators within the frozen numerical policy."]; cert["claim"]["prohibited"]=["Claim that prediction is impossible under every restricted population model."]
    else:
        cert["claim"]["permitted"]=["Numerical status unresolved; abstain."]; cert["claim"]["prohibited"]=["ATI0 exact claim.","ATI5 no-go claim without stronger numerical evidence."]
    return cert
