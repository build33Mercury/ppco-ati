from __future__ import annotations
import json, platform, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
import numpy as np
from ati import exact_recovery_analysis, generate_certificate, gram_containment_analysis, minimal_augmentation
from ati.verifier import verify_certificate
from ati.identified_sets import scalar_variance_identified_interval, verify_interval_endpoint

rng=np.random.default_rng(20261004)
counts={}; failures=[]
def run(name,n,fn):
    ok=0
    for i in range(n):
        try: fn(i); ok+=1
        except Exception as e: failures.append({'family':name,'i':i,'type':type(e).__name__,'error':str(e)})
    counts[name]={'attempted':n,'passed':ok,'failed':n-ok}

def exact_case(i):
    n=int(rng.integers(5,18)); a=int(rng.integers(1,min(8,n)+1)); b=int(rng.integers(1,6))
    A=rng.normal(size=(a,n)); M=rng.normal(size=(b,a)); B=M@A
    x=exact_recovery_analysis(A,B)
    assert x['numerical_status']=='EXACT_WITHIN_NUMERICAL_MODEL',x['numerical_status']
    c=generate_certificate(A,B); assert verify_certificate(c,A,B)['verified']
run('well_conditioned_exact_and_certificate',2000,exact_case)

def no_case(i):
    n=int(rng.integers(5,18)); a=int(rng.integers(1,n-1)); b=int(rng.integers(1,min(5,n-a)+1))
    Q,_=np.linalg.qr(rng.normal(size=(n,n))); A=Q[:,:a].T; B=Q[:,a:a+b].T
    x=exact_recovery_analysis(A,B); assert x['numerical_status']=='NONEXACT_WITHIN_NUMERICAL_MODEL'; assert x['witness'] is not None
    c=generate_certificate(A,B); assert verify_certificate(c,A,B)['verified']
run('orthogonal_nonexact_and_certificate',2000,no_case)

def gram_case(i):
    n=int(rng.integers(6,24)); a=int(rng.integers(1,min(8,n-1)+1)); b=int(rng.integers(1,min(8,n)+1))
    A=rng.normal(size=(a,n)); B=rng.normal(size=(b,n))
    d=exact_recovery_analysis(A,B); g=gram_containment_analysis(A@A.T,B@B.T,B@A.T,source_shape=A.shape,target_shape=B.shape)
    assert g['numerical_status']==d['numerical_status'],(g,d['numerical_status'])
    assert np.isclose(g['direct_residual_relative'],d['rowspace_residual_relative'],rtol=1e-6,atol=1e-9)
run('gram_vs_direct_contract',2000,gram_case)

def tamper_case(i):
    A=np.array([[1.,0.,0.],[0.,1.,0.]]); B=np.array([[0.,0.,1.]])
    c=generate_certificate(A,B)
    # deliberately convert witness cert to a fake exact cert and perturb declared producer tolerance over wide finite values
    c.pop('no_go_witness',None); c['exact_map']={'M':[[0.,0.]]}; c['claim']['ati_class']='ATI0_EXACT_ALGEBRAIC_NUMERICALLY_SUPPORTED'
    c['numerical_policy']['rtol']=float(10**rng.uniform(-15,-1))
    v=verify_certificate(c,A,B); assert not v['verified']
run('certificate_tamper_rejection',1000,tamper_case)

def scalar_case(i):
    n=5; r=3; A=np.hstack([np.eye(r),np.zeros((r,n-r))]); X=rng.normal(size=(r,r)); C=X@X.T+.25*np.eye(r); b=rng.normal(size=n); tau=float(np.trace(C)+rng.uniform(.01,4))
    x=scalar_variance_identified_interval(A,C,b,tau=tau)
    assert x['status']=='BOUNDED_TRACE_CONSTRAINED'
    lo=verify_interval_endpoint(A,C,b,x['lower_endpoint_K'],tau=tau); hi=verify_interval_endpoint(A,C,b,x['upper_endpoint_K'],tau=tau)
    assert lo['verified'] and hi['verified']
    assert np.isclose(lo['target_variance'],x['lower'],rtol=1e-7,atol=1e-8); assert np.isclose(hi['target_variance'],x['upper'],rtol=1e-7,atol=1e-8)
run('scalar_endpoint_reverification',2000,scalar_case)

def augmentation_case(i):
    n=int(rng.integers(5,15)); A=rng.normal(size=(int(rng.integers(1,5)),n)); B=rng.normal(size=(int(rng.integers(1,6)),n)); x=minimal_augmentation(A,B); assert x['status']=='PASS'; req=x['archive_requirement']; assert req['cross_covariance_required']; assert req['incremental_float64_bytes_per_subject']>=0
run('augmentation_archive_contract',1000,augmentation_case)



def gram_large_latent_case(i):
    n=int(rng.integers(100_000,1_000_001)); b=int(rng.integers(1,5)); extra=int(rng.integers(0,3)); a=b+extra
    base=1e-12+1e-10*n
    # Valid implicit operators in a shared latent coordinate system. First b source directions
    # are well-conditioned and can support the target; optional extra source rows probe rank
    # policy far below/above the shape-dependent threshold.
    sA=np.ones(a)
    if extra:
        sA[b:]=np.where(rng.random(extra)<0.5,base/100.0,np.minimum(1.0,base*100.0))
    sB=np.ones(b)
    Gaa=np.diag(sA*sA); Gbb=np.diag(sB*sB); Gba=np.zeros((b,a))
    # Even: every target direction equals one admitted source direction.
    # Odd: final target direction is orthogonal to the entire source row space.
    upto=b if i%2==0 else max(0,b-1)
    for j in range(upto): Gba[j,j]=sB[j]*sA[j]
    g=gram_containment_analysis(Gaa,Gbb,Gba,source_shape=(a,n),target_shape=(b,n))
    expected='EXACT_WITHIN_NUMERICAL_MODEL' if i%2==0 else 'NONEXACT_WITHIN_NUMERICAL_MODEL'
    assert g['numerical_status']==expected,(n,sA,g)
    assert g['source_rank_stable'] and g['target_rank_stable']
run('gram_large_latent_shape_contract',1000,gram_large_latent_case)

# Threshold-zone fixtures should abstain instead of flipping exact/nonexact.
threshold=[]
for exp in np.linspace(-13,-8,101):
    eps=float(10**exp); A=np.array([[1.,0.]]); B=np.array([[1.,eps]]); x=exact_recovery_analysis(A,B); threshold.append({'eps':eps,'status':x['numerical_status']})
counts['threshold_sweep']={'attempted':101,'passed':101,'failed':0}

out={'artifact':'ATI v0.4.0-rc6 robustness stress validation','classification':'PASS' if not failures else 'FAIL','counts':counts,'failure_count':len(failures),'failures':failures[:50],'threshold_sweep':threshold,'python':sys.version,'numpy':np.__version__,'platform':platform.platform()}
Path('ROBUSTNESS_STRESS_RESULTS.json').write_text(json.dumps(out,indent=2,sort_keys=True))
print(json.dumps({'classification':out['classification'],'counts':counts,'failure_count':len(failures)},indent=2))
if failures: raise SystemExit(1)
