from __future__ import annotations
from typing import Any, Dict, Optional
import numpy as np

from .core import _as_matrix
from .numerics import (
    DEFAULT_RTOL, DEFAULT_ATOL, DEFAULT_STABILITY_FACTOR,
    validate_tolerances, absrel_threshold, tight_loose, svd_rank_info, pinv_from_svd, policy_record,
)


def _checked_symmetric(C,name,rtol,atol,stability_factor):
    C=np.asarray(C,dtype=np.float64)
    if C.ndim!=2 or C.shape[0]!=C.shape[1] or not np.all(np.isfinite(C)): raise ValueError(f"{name} must be a finite square matrix")
    asym=float(np.linalg.norm(C-C.T,2)) if C.size else 0.0; scale=float(np.linalg.norm(C,2)) if C.size else 0.0
    base=absrel_threshold(scale,rtol=rtol,atol=atol); tight,loose=tight_loose(base,stability_factor)
    if asym>loose: raise ValueError(f"{name} is materially asymmetric")
    status="SYMMETRIC_WITHIN_NUMERICAL_MODEL" if asym<=tight else "UNKNOWN_NUMERICAL"
    return 0.5*(C+C.T),status,asym,tight,loose


def _orthogonal_completion(Qs,n):
    if Qs.shape[0]==0: return np.eye(n)
    _,_,vh=np.linalg.svd(Qs,full_matrices=True); return vh[Qs.shape[0]:,:]


def _source_block(A,C_A,rtol,atol,stability_factor):
    A=_as_matrix(A,"A"); C_A,sym_status,asym,sym_tight,sym_loose=_checked_symmetric(C_A,"C_A",rtol,atol,stability_factor)
    if C_A.shape!=(A.shape[0],A.shape[0]): raise ValueError("C_A shape must match source outputs")
    ia=svd_rank_info(A,rtol=rtol,atol=atol,stability_factor=stability_factor,full_matrices=False)
    if not ia["rank_stable"]:
        return {"A":A,"C_A":C_A,"source_status":"UNKNOWN_NUMERICAL","source_feasible":None,"reason":"A_rank_unstable","symmetry_status":sym_status,"asymmetry_norm2":asym}
    Qs=ia["vh"][:ia["rank"],:]; Qh=_orthogonal_completion(Qs,A.shape[1]); L=A@Qs.T
    if ia["rank"]==0:
        Kss=np.zeros((0,0)); recon=np.zeros_like(C_A); il={"rank_stable":True}
    else:
        Lp,il=pinv_from_svd(L,rtol=rtol,atol=atol,stability_factor=stability_factor,require_stable_rank=True)
        if Lp is None: return {"A":A,"C_A":C_A,"source_status":"UNKNOWN_NUMERICAL","source_feasible":None,"reason":"L_rank_unstable","symmetry_status":sym_status,"asymmetry_norm2":asym}
        Kss=0.5*((Lp@C_A@Lp.T)+(Lp@C_A@Lp.T).T); recon=L@Kss@L.T
    scale=max(float(np.linalg.norm(C_A,2)) if C_A.size else 0.0,1.0); base=absrel_threshold(scale,rtol=rtol,atol=atol); tight,loose=tight_loose(base,stability_factor)
    resid=float(np.linalg.norm(recon-C_A,2)) if C_A.size else 0.0; eig=np.linalg.eigvalsh(Kss) if Kss.size else np.array([]); mine=float(eig.min()) if eig.size else 0.0
    if resid>loose or mine < -loose:
        status="INFEASIBLE_SOURCE_COVARIANCE"; feasible=False
    elif sym_status=="UNKNOWN_NUMERICAL" or (tight<resid<=loose) or (-loose<=mine < -tight):
        status="UNKNOWN_NUMERICAL"; feasible=None
    else:
        status="FEASIBLE"; feasible=True
        if eig.size and mine<0:
            vals,vecs=np.linalg.eigh(Kss); vals[np.abs(vals)<=tight]=0.0; Kss=(vecs*vals)@vecs.T; Kss=0.5*(Kss+Kss.T)
    return {"A":A,"C_A":C_A,"Qs":Qs,"Qh":Qh,"L":L,"Kss":Kss,"source_status":status,"source_feasible":feasible,
            "source_consistency_residual_norm2":resid,"source_psd_min_eigenvalue":mine,"source_threshold_tight":tight,"source_threshold_loose":loose,
            "symmetry_status":sym_status,"asymmetry_norm2":asym}


def _construct_scalar_endpoint(block,b,hidden_amplitude,sign):
    Qs,Qh,Kss=block["Qs"],block["Qh"],block["Kss"]; n=block["A"].shape[1]; b=np.asarray(b,dtype=np.float64).reshape(-1)
    bs=Qs@b if Qs.size else np.zeros(0); bh=Qh@b if Qh.size else np.zeros(0); hnorm=float(np.linalg.norm(bh)); a2=float(bs@Kss@bs) if bs.size else 0.0; a=float(np.sqrt(max(a2,0.0))); s=float(hidden_amplitude)
    Ksh=np.zeros((Kss.shape[0],Qh.shape[0])); Khh=np.zeros((Qh.shape[0],Qh.shape[0]))
    if hnorm>0 and s>0:
        u=bh/hnorm; alpha=(s/hnorm)**2; Khh=alpha*np.outer(u,u)
        if a>0 and bs.size:
            g=(Kss@bs)/a; Ksh=sign*(s/hnorm)*np.outer(g,u)
    top=np.concatenate([Kss,Ksh],axis=1) if Qh.shape[0] else Kss
    if Qh.shape[0]:
        bottom=np.concatenate([Ksh.T,Khh],axis=1); Kbasis=np.concatenate([top,bottom],axis=0); O=np.concatenate([Qs,Qh],axis=0)
    else: Kbasis=Kss; O=Qs
    K=O.T@Kbasis@O if O.size else np.zeros((n,n)); return 0.5*(K+K.T)


def scalar_variance_identified_interval(A,C_A,b,tau:Optional[float]=None,rtol=DEFAULT_RTOL,atol=DEFAULT_ATOL,
                                        stability_factor=DEFAULT_STABILITY_FACTOR,return_endpoint_witnesses=True)->Dict[str,Any]:
    rtol,atol,stability_factor=validate_tolerances(rtol,atol,stability_factor)
    # Validate tau BEFORE any exact shortcut.
    if tau is not None:
        tau=float(tau)
        if not np.isfinite(tau) or tau<0: raise ValueError("tau must be finite and nonnegative")
    block=_source_block(A,C_A,rtol,atol,stability_factor); A=block["A"]; b=np.asarray(b,dtype=np.float64).reshape(-1)
    if b.shape[0]!=A.shape[1] or not np.all(np.isfinite(b)): raise ValueError("b must be finite and match latent dimension")
    base={"estimand":"scalar_target_variance_b_K_bT","source_feasible":block.get("source_feasible"),"source_status":block.get("source_status"),
          "source_consistency_residual_norm2":block.get("source_consistency_residual_norm2"),"source_psd_min_eigenvalue":block.get("source_psd_min_eigenvalue"),
          "source_covariance_symmetry_status":block.get("symmetry_status"),"source_covariance_asymmetry_norm2":block.get("asymmetry_norm2"),
          "tau":tau,"numerical_policy":policy_record(rtol=rtol,atol=atol,stability_factor=stability_factor),"rtol":rtol,"atol":atol,"stability_factor":stability_factor}
    if block.get("source_status")=="INFEASIBLE_SOURCE_COVARIANCE": return {**base,"status":"INFEASIBLE_SOURCE_COVARIANCE","lower":None,"upper":None,"upper_unbounded":False}
    if block.get("source_status")!="FEASIBLE": return {**base,"status":"UNKNOWN_NUMERICAL","lower":None,"upper":None,"upper_unbounded":False,"reason":block.get("reason","source_feasibility_ambiguous")}
    Qs,Qh,Kss=block["Qs"],block["Qh"],block["Kss"]; bs=Qs@b if Qs.size else np.zeros(0); bh=Qh@b if Qh.size else np.zeros(0)
    a2=max(float(bs@Kss@bs) if bs.size else 0.0,0.0); a=float(np.sqrt(a2)); hnorm=float(np.linalg.norm(bh)); tr=float(np.trace(Kss))
    base.update({"identified_source_component_variance":a2,"target_hidden_component_norm2":hnorm,"minimum_latent_trace_consistent_with_source":tr})
    bbase=absrel_threshold(float(np.linalg.norm(b)),rtol=rtol,atol=atol); _,hloose=tight_loose(bbase,stability_factor)
    # A target-null component is substantively different from an ordinary residual: do not collapse a small but resolved hidden direction merely because rtol is loose. Exactness requires it to be below the absolute numerical floor; values above that floor but below the robust nonzero threshold are UNKNOWN.
    htight=max(float(atol), 64.0*np.finfo(float).eps*max(float(np.linalg.norm(b)),1.0))
    hidden_status="ZERO_WITHIN_NUMERICAL_MODEL" if hnorm<=htight else "NONZERO_WITHIN_NUMERICAL_MODEL" if hnorm>hloose else "UNKNOWN_NUMERICAL"
    base["hidden_component_status"]=hidden_status
    if hidden_status=="UNKNOWN_NUMERICAL": return {**base,"status":"UNKNOWN_NUMERICAL","lower":None,"upper":None,"upper_unbounded":False,"reason":"hidden_component_near_threshold"}

    slack=None
    if tau is not None:
        slack=float(tau-tr); tbase=absrel_threshold(max(tau,tr),rtol=rtol,atol=atol); ttight,tloose=tight_loose(tbase,stability_factor)
        if slack < -tloose: return {**base,"status":"INFEASIBLE_TRACE_BOUND","lower":None,"upper":None,"upper_unbounded":False,"remaining_trace_budget":slack}
        if -tloose <= slack < -ttight: return {**base,"status":"UNKNOWN_NUMERICAL","lower":None,"upper":None,"upper_unbounded":False,"remaining_trace_budget":slack,"reason":"trace_feasibility_near_threshold"}
        slack=max(slack,0.0)

    if hidden_status=="ZERO_WITHIN_NUMERICAL_MODEL":
        K=_construct_scalar_endpoint(block,b,0,+1); result={**base,"status":"EXACT_WITHIN_NUMERICAL_MODEL","lower":a2,"upper":a2,"upper_unbounded":False,"remaining_trace_budget":slack}
        if return_endpoint_witnesses: result["lower_endpoint_K"]=K; result["upper_endpoint_K"]=K
        return result
    if tau is None:
        Klo=_construct_scalar_endpoint(block,b,a,-1); result={**base,"status":"UNBOUNDED_UPPER","lower":0.0,"upper":None,"upper_unbounded":True,"remaining_trace_budget":None,
                                                                  "warning":"No finite upper bound under unrestricted PSD K when target has a source-null component."}
        if return_endpoint_witnesses: result["lower_endpoint_K"]=Klo
        return result
    smax=float(np.sqrt(slack)*hnorm); lower_s=min(a,smax); lower=float((a-lower_s)**2); upper=float((a+smax)**2)
    Klo=_construct_scalar_endpoint(block,b,lower_s,-1); Khi=_construct_scalar_endpoint(block,b,smax,+1)
    result={**base,"status":"BOUNDED_TRACE_CONSTRAINED","lower":lower,"upper":upper,"upper_unbounded":False,"remaining_trace_budget":slack,
            "hidden_standard_deviation_budget":smax,"normalized_width_over_max_abs_endpoint":float((upper-lower)/max(abs(upper),abs(lower),atol)),
            "warning":"Conditional on the declared trace(K)<=tau constraint; identified intervals/PSD-fiber optimization have prior art, so no generic novelty claim is permitted."}
    if return_endpoint_witnesses: result["lower_endpoint_K"]=Klo; result["upper_endpoint_K"]=Khi
    return result


def identified_interval_record(result):
    out={k:v for k,v in result.items() if not k.endswith("_K")}; out["schema_version"]="ATI-scalar-interval-v0.2"
    for k,v in list(out.items()):
        if isinstance(v,np.generic): out[k]=v.item()
    return out


def verify_interval_endpoint(A,C_A,b,K,tau=None,rtol=DEFAULT_RTOL,atol=DEFAULT_ATOL,stability_factor=DEFAULT_STABILITY_FACTOR):
    rtol,atol,stability_factor=validate_tolerances(rtol,atol,stability_factor,verifier=True)
    if tau is not None:
        tau=float(tau)
        if not np.isfinite(tau) or tau<0: raise ValueError("tau must be finite and nonnegative")
    A=_as_matrix(A,"A"); C_A,cs,ca,_,_=_checked_symmetric(C_A,"C_A",rtol,atol,stability_factor); b=np.asarray(b,dtype=np.float64).reshape(-1); K=np.asarray(K,dtype=np.float64)
    if K.shape!=(A.shape[1],A.shape[1]) or not np.all(np.isfinite(K)): raise ValueError("K has wrong shape or non-finite values")
    Ks,ks,ka,_,_=_checked_symmetric(K,"K",rtol,atol,stability_factor)
    scale=max(float(np.linalg.norm(C_A,2)) if C_A.size else 0.0,1.0); base=absrel_threshold(scale,rtol=rtol,atol=atol); tight,loose=tight_loose(base,stability_factor)
    eig=np.linalg.eigvalsh(Ks); mine=float(eig.min()) if eig.size else 0.0; srcres=float(np.linalg.norm(A@Ks@A.T-C_A,2)); tr=float(np.trace(Ks)); val=float(b@Ks@b)
    issues=[]; unknown=[]
    if ks=="UNKNOWN_NUMERICAL" or cs=="UNKNOWN_NUMERICAL": unknown.append("symmetry_near_threshold")
    if mine < -loose: issues.append("K_not_PSD")
    elif mine < -tight: unknown.append("K_PSD_near_threshold")
    if srcres>loose: issues.append("source_constraint_failed")
    elif srcres>tight: unknown.append("source_constraint_near_threshold")
    if tau is not None:
        tbase=absrel_threshold(max(tau,1.0),rtol=rtol,atol=atol); ttight,tloose=tight_loose(tbase,stability_factor); excess=tr-tau
        if excess>tloose: issues.append("trace_bound_failed")
        elif excess>ttight: unknown.append("trace_bound_near_threshold")
    status="REJECTED" if issues else "UNKNOWN_NUMERICAL" if unknown else "VERIFIED"
    return {"verified":status=="VERIFIED","status":status,"issues":issues,"unknown_reasons":unknown,"target_variance":val,"trace":tr,"min_eigenvalue":mine,"source_residual_norm2":srcres,"K_asymmetry_norm2":ka,"C_A_asymmetry_norm2":ca}
