from __future__ import annotations
import json, platform, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
import numpy as np

from ati import exact_recovery_analysis, minimal_augmentation, universal_representation, trace_bounded_projection_error_bound
from ati.identified_sets import scalar_variance_identified_interval, verify_interval_endpoint
from ati.exact_verifier import exact_rational_classify
from ati.stability import perturbation_stability_sweep

RNG=np.random.default_rng(20261001)
counts={}
failures=[]

def run(name,n,fn):
    ok=0
    for i in range(n):
        try:
            fn(i); ok+=1
        except Exception as e:
            failures.append({'family':name,'index':i,'type':type(e).__name__,'error':str(e)})
    counts[name]={'passed':ok,'attempted':n,'failed':n-ok}

run('exact_nested',500,lambda i: (
    (lambda A,M,B: (_ for _ in ()).throw(AssertionError('not exact')) if exact_recovery_analysis(A,B)['numerical_status']!='EXACT_WITHIN_NUMERICAL_MODEL' else None)(
        (A:=RNG.normal(size=(4,9))), (M:=RNG.normal(size=(3,4))), M@A)
))

def nonnested(_):
    Q,_=np.linalg.qr(RNG.normal(size=(9,9)))
    A=Q[:,:3].T; B=Q[:,3:5].T
    r=exact_recovery_analysis(A,B)
    assert r['numerical_status']=='NONEXACT_WITHIN_NUMERICAL_MODEL' and r['witness'] is not None
run('nonnested',500,nonnested)

def aug(_):
    A=RNG.normal(size=(3,9)); B=RNG.normal(size=(4,9))
    x=minimal_augmentation(A,B)
    assert x['verification_status']=='EXACT_WITHIN_NUMERICAL_MODEL'
run('minimal_augmentation',500,aug)

def uni(_):
    ops=[RNG.normal(size=(2,8)),RNG.normal(size=(3,8)),RNG.normal(size=(1,8))]
    x=universal_representation(ops)
    assert all(s=='EXACT_WITHIN_NUMERICAL_MODEL' for s in x['target_verification_statuses'])
run('universal_representation',200,uni)

def trace_bound(_):
    n=8; A=RNG.normal(size=(3,n)); B=RNG.normal(size=(4,n)); X=RNG.normal(size=(n,n)); K=X@X.T
    out=trace_bounded_projection_error_bound(A,B,float(np.trace(K)))
    _,_,vh=np.linalg.svd(A,full_matrices=False); rank=np.linalg.matrix_rank(A); P=vh[:rank].T@vh[:rank]
    lhs=float(np.linalg.norm(B@K@B.T-B@P@K@P@B.T,2))
    assert lhs <= out['bound'] + 1e-8*max(1.0,out['bound'])
run('trace_projection_bound',2000,trace_bound)

def interval_completion(_):
    # Coordinate source so feasible completions can be sampled exactly.
    n=5; r=3
    A=np.hstack([np.eye(r),np.zeros((r,n-r))])
    X=RNG.normal(size=(r,r)); C=X@X.T + .2*np.eye(r)
    b=RNG.normal(size=n)
    tau=float(np.trace(C)+RNG.uniform(.1,4.0))
    res=scalar_variance_identified_interval(A,C,b,tau=tau)
    assert res['status']=='BOUNDED_TRACE_CONSTRAINED'
    lo=verify_interval_endpoint(A,C,b,res['lower_endpoint_K'],tau=tau)
    hi=verify_interval_endpoint(A,C,b,res['upper_endpoint_K'],tau=tau)
    assert lo['verified'] and hi['verified']
    assert abs(lo['target_variance']-res['lower']) <= 1e-7*max(1.0,abs(res['lower']))
    assert abs(hi['target_variance']-res['upper']) <= 1e-7*max(1.0,abs(res['upper']))
    # Random PSD completion with same source block via factor stacking.
    vals,vecs=np.linalg.eigh(C); Ls=vecs@np.diag(np.sqrt(np.maximum(vals,0)))
    slack=tau-np.trace(C)
    H=RNG.normal(size=(n-r,r))
    norm2=float(np.sum(H*H))
    if norm2:
        H=H*np.sqrt(slack*RNG.uniform(0,1)/norm2)
    F=np.vstack([Ls,H]); K=F@F.T
    val=float(b@K@b)
    assert res['lower']-1e-7 <= val <= res['upper']+1e-7
run('scalar_identified_interval',1000,interval_completion)

def exact_rat(_):
    # Random integer exact containment.
    A=RNG.integers(-3,4,size=(3,5)).tolist()
    M=RNG.integers(-2,3,size=(2,3))
    B=(M@np.asarray(A,dtype=int)).tolist()
    x=exact_rational_classify(A,B)
    assert x['rowspace_contained'] and x['verification']=='B_EQUALS_MA_EXACTLY'
run('exact_rational_containment',200,exact_rat)

# deterministic perturbation smoke: exact baseline and nonexact baseline
try:
    exact_sweep=perturbation_stability_sweep(np.array([[1.,0.,0.],[0.,1.,0.]]),np.array([[1.,1.,0.]]),epsilons=[0,1e-12,1e-10,1e-8],trials=25,seed=20261001)
    nonexact_sweep=perturbation_stability_sweep(np.array([[1.,0.,0.]]),np.array([[0.,1.,0.]]),epsilons=[0,1e-12,1e-10,1e-8],trials=25,seed=20261002,perturb_source=False)
    counts['perturbation_sweeps']={'passed':2,'attempted':2,'failed':0}
except Exception as e:
    exact_sweep=None; nonexact_sweep=None
    counts['perturbation_sweeps']={'passed':0,'attempted':2,'failed':2}
    failures.append({'family':'perturbation_sweeps','type':type(e).__name__,'error':str(e)})

result={
 'artifact':'ATI methods branch v0.4.0-rc6 adversarial development audit',
 'date':'2026-10-04',
 'classification':'PASS' if not failures else 'FAIL',
 'development_only':True,
 'counts':counts,
 'failures':failures,
 'perturbation_exact_fixture':exact_sweep,
 'perturbation_nonexact_fixture':nonexact_sweep,
 'python':sys.version,
 'platform':platform.platform(),
 'numpy':np.__version__,
 'interpretation':'Randomized numerical development evidence only. Not formal proof, real-atlas validation, or manuscript outcome authority.'
}
Path('DEVELOPMENT_AUDIT_V04RC6.json').write_text(json.dumps(result,indent=2,sort_keys=True),encoding='utf-8')
print(json.dumps({'classification':result['classification'],'counts':counts,'failure_count':len(failures)},indent=2))
if failures: raise SystemExit(1)
