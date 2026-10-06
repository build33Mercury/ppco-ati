from __future__ import annotations
from typing import Any, Dict, Iterable
import numpy as np
from .core import _as_matrix, exact_recovery_analysis
from .numerics import DEFAULT_RTOL, DEFAULT_ATOL, DEFAULT_STABILITY_FACTOR

STATUSES=['EXACT_WITHIN_NUMERICAL_MODEL','NONEXACT_WITHIN_NUMERICAL_MODEL','UNKNOWN_NUMERICAL']

def _fro_normalized_noise(rng,shape):
    E=rng.normal(size=shape); n=float(np.linalg.norm(E,'fro')); return E/n if n else E

def perturbation_stability_sweep(A,B,epsilons:Iterable[float]=(0.0,1e-12,1e-10,1e-8,1e-6),trials:int=20,seed:int=20261001,
                                 perturb_source:bool=True,perturb_target:bool=True,rtol:float=DEFAULT_RTOL,atol:float=DEFAULT_ATOL,
                                 stability_factor:float=DEFAULT_STABILITY_FACTOR)->Dict[str,Any]:
    A=_as_matrix(A,'A'); B=_as_matrix(B,'B')
    if A.shape[1]!=B.shape[1]: raise ValueError('A and B must act on same latent dimension')
    if trials<1: raise ValueError('trials must be >=1')
    eps=[float(x) for x in epsilons]
    if any((not np.isfinite(x) or x<0) for x in eps): raise ValueError('epsilons must be finite and nonnegative')
    rng=np.random.default_rng(seed); normA=max(float(np.linalg.norm(A,'fro')),1.0); normB=max(float(np.linalg.norm(B,'fro')),1.0); rows=[]
    for e in eps:
        counts={k:0 for k in STATUSES}; residuals=[]; ntrials=1 if e==0 else trials
        for _ in range(ntrials):
            Ap=A.copy(); Bp=B.copy()
            if e>0 and perturb_source: Ap=Ap+e*normA*_fro_normalized_noise(rng,A.shape)
            if e>0 and perturb_target: Bp=Bp+e*normB*_fro_normalized_noise(rng,B.shape)
            x=exact_recovery_analysis(Ap,Bp,rtol=rtol,atol=atol,stability_factor=stability_factor); counts[x['numerical_status']]+=1
            if x['rowspace_residual_relative'] is not None: residuals.append(x['rowspace_residual_relative'])
        rows.append({'epsilon_relative_fro':e,'trials':ntrials,'status_counts':counts,
                     'residual_relative_min':float(min(residuals)) if residuals else None,
                     'residual_relative_median':float(np.median(residuals)) if residuals else None,
                     'residual_relative_max':float(max(residuals)) if residuals else None})
    observed={k for row in rows for k,v in row['status_counts'].items() if v}
    return {'seed':int(seed),'perturb_source':bool(perturb_source),'perturb_target':bool(perturb_target),'rtol':float(rtol),'atol':float(atol),'stability_factor':float(stability_factor),'rows':rows,'observed_statuses':sorted(observed),'warning':'Generic numerical/operator perturbation stress test; not calibrated registration uncertainty.'}
