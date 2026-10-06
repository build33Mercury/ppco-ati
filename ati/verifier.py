from __future__ import annotations

from typing import Dict, Any
import hashlib
import numpy as np

from .numerics import (
    DEFAULT_RTOL, DEFAULT_ATOL, DEFAULT_STABILITY_FACTOR,
    validate_tolerances, absrel_threshold, tight_loose, svd_rank_info, pinv_from_svd,
)


def _mat(X):
    X=np.asarray(X,dtype=np.float64)
    if X.ndim!=2 or not np.all(np.isfinite(X)): raise ValueError("invalid matrix")
    return np.ascontiguousarray(X)


def _hash(X):
    X=_mat(X); h=hashlib.sha256(); h.update(str(X.shape).encode("ascii")); h.update(b"|float64|C|"); h.update(X.astype("<f8",copy=False).tobytes(order="C")); return h.hexdigest()


def _local_status(A,B,rtol,atol,stability_factor):
    ia=svd_rank_info(A,rtol=rtol,atol=atol,stability_factor=stability_factor,full_matrices=False)
    if not ia["rank_stable"]: return "UNKNOWN_NUMERICAL", None, None
    Q=ia["vh"][:ia["rank"],:]; P=np.zeros((A.shape[1],A.shape[1])) if ia["rank"]==0 else Q.T@Q
    R=B@(np.eye(B.shape[1])-P); nr=float(np.linalg.norm(R,2)) if R.size else 0.0; nb=float(np.linalg.norm(B,2)) if B.size else 0.0
    base=absrel_threshold(nb,rtol=rtol,atol=atol); tight,loose=tight_loose(base,stability_factor)
    if nr<=tight: return "EXACT_WITHIN_NUMERICAL_MODEL",nr,(tight,loose)
    if nr>loose: return "NONEXACT_WITHIN_NUMERICAL_MODEL",nr,(tight,loose)
    return "UNKNOWN_NUMERICAL",nr,(tight,loose)


def verify_certificate(cert: Dict[str,Any], A, B, *,
                       rtol: float=DEFAULT_RTOL, atol: float=DEFAULT_ATOL,
                       stability_factor: float=DEFAULT_STABILITY_FACTOR) -> Dict[str,Any]:
    """Independent verifier. Producer-provided tolerances are diagnostics only and never control acceptance."""
    rtol,atol,stability_factor=validate_tolerances(rtol,atol,stability_factor,verifier=True)
    A=_mat(A); B=_mat(B); issues=[]
    if A.shape[1]!=B.shape[1]: raise ValueError("A and B latent dimensions differ")
    schema=cert.get("schema_version")
    if schema!="ATI-certificate-v0.3": issues.append("unsupported_schema_version")
    if cert.get("source",{}).get("operator_sha256")!=_hash(A): issues.append("source_hash_mismatch")
    if cert.get("target",{}).get("operator_sha256")!=_hash(B): issues.append("target_hash_mismatch")
    if cert.get("source",{}).get("shape")!=list(A.shape): issues.append("source_shape_mismatch")
    if cert.get("target",{}).get("shape")!=list(B.shape): issues.append("target_shape_mismatch")
    if "exact_map" in cert and "no_go_witness" in cert: issues.append("conflicting_payloads")

    local_status, local_resid, thresholds=_local_status(A,B,rtol,atol,stability_factor)
    declared=cert.get("claim",{}).get("ati_class")
    outbase={"local_numerical_status":local_status,"verifier_rtol":rtol,"verifier_atol":atol,"verifier_stability_factor":stability_factor,"producer_tolerances_ignored_for_acceptance":True}

    if "exact_map" in cert:
        M=np.asarray(cert["exact_map"].get("M"),dtype=np.float64)
        if M.ndim!=2 or M.shape!=(B.shape[0],A.shape[0]) or not np.all(np.isfinite(M)):
            issues.append("exact_map_shape_or_values_invalid"); return {**outbase,"verified":False,"mode":"exact_map","issues":issues}
        resid=float(np.linalg.norm(B-M@A,2)); nb=float(np.linalg.norm(B,2)) if B.size else 0.0
        _, loose=tight_loose(absrel_threshold(nb,rtol=rtol,atol=atol),stability_factor)
        if local_status!="EXACT_WITHIN_NUMERICAL_MODEL": issues.append("independent_rowspace_test_not_exact")
        if resid>loose: issues.append("exact_map_residual_too_large")
        if declared and not declared.startswith("ATI0_"): issues.append("declared_class_inconsistent_with_exact_payload")
        return {**outbase,"verified":len(issues)==0,"mode":"exact_map","residual_norm2":resid,"acceptance_threshold_loose":loose,"issues":issues}

    if "no_go_witness" in cert:
        v=np.asarray(cert["no_go_witness"].get("v"),dtype=np.float64)
        if v.ndim!=1 or v.shape[0]!=A.shape[1] or not np.all(np.isfinite(v)):
            issues.append("witness_shape_or_values_invalid"); return {**outbase,"verified":False,"mode":"witness","issues":issues}
        av=float(np.linalg.norm(A@v)); bv=float(np.linalg.norm(B@v)); nv=float(np.linalg.norm(v))
        _, looseA=tight_loose(absrel_threshold(float(np.linalg.norm(A,2)) if A.size else 0.0,rtol=rtol,atol=atol),stability_factor)
        _, looseB=tight_loose(absrel_threshold(float(np.linalg.norm(B,2)) if B.size else 0.0,rtol=rtol,atol=atol),stability_factor)
        if local_status!="NONEXACT_WITHIN_NUMERICAL_MODEL": issues.append("independent_rowspace_test_not_nonexact")
        if av>looseA: issues.append("witness_not_in_kernel")
        if bv<=looseB: issues.append("witness_not_visible_to_target")
        if abs(nv-1.0)>max(1e-10,10*atol): issues.append("witness_not_unit_norm")
        if declared and not declared.startswith("ATI5_"): issues.append("declared_class_inconsistent_with_witness_payload")
        return {**outbase,"verified":len(issues)==0,"mode":"witness","A_v_norm2":av,"B_v_norm2":bv,"v_norm2":nv,"issues":issues}

    producer_status=cert.get("numerics",{}).get("status")
    if producer_status!="UNKNOWN_NUMERICAL": issues.append("no_verifiable_payload")
    if declared not in {None,"UNKNOWN_NUMERICAL"}: issues.append("declared_class_inconsistent_with_no_payload")
    # Conservative abstention can be accepted even if this verifier could decide more strongly.
    return {**outbase,"verified":len(issues)==0,"mode":"unknown_abstention","issues":issues}
