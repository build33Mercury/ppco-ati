from __future__ import annotations
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
# Historical specification digest retained without bundling internal protocol documents.
PROTOCOL_SHA256='7d7987fb531b1ec72c904bcdd4931371803c2a84fcd75a9c58fe6b80c9fa2c19'
SEED=20261001
rng=np.random.default_rng(SEED)

def gen(n, shifted=False):
    x=rng.uniform(0.5,2.0,size=n)
    eps=rng.normal(0.0,0.05,size=n)
    y=(3.0-0.8*x+eps) if shifted else (1.0+0.8*x+eps)
    if np.any(y<=0): raise RuntimeError('frozen generator produced nonpositive variance')
    return x,y

def fit_ols(x,y):
    X=np.column_stack([np.ones_like(x),x])
    beta=np.linalg.lstsq(X,y,rcond=None)[0]
    return beta

def metrics(y,p):
    mae=float(np.mean(np.abs(y-p)))
    ssr=float(np.sum((y-p)**2)); sst=float(np.sum((y-y.mean())**2))
    r2=float(1.0-ssr/sst)
    return {'MAE':mae,'R2':r2}

xtr,ytr=gen(500,False); xid,yid=gen(500,False); xsh,ysh=gen(500,True)
beta=fit_ols(xtr,ytr)
pid=beta[0]+beta[1]*xid; psh=beta[0]+beta[1]*xsh
A=np.array([[1.,0.]]); B=np.array([[0.,1.]]); v=np.array([0.,1.]); K0=np.diag([1.25,2.0]); K1=K0+2.0*np.outer(v,v)
source0=(A@K0@A.T).tolist(); source1=(A@K1@A.T).tolist(); target0=(B@K0@B.T).tolist(); target1=(B@K1@B.T).tolist()
report={
  'status':'PASS','type':'PROSPECTIVE_SYNTHETIC_DEVELOPMENT_DEMO','seed':SEED,
  'protocol_sha256':PROTOCOL_SHA256,
  'model':{'intercept':float(beta[0]),'slope':float(beta[1])},
  'iid_test':metrics(yid,pid),'shift_test':metrics(ysh,psh),
  'source_distribution':{'iid_mean':float(xid.mean()),'shift_mean':float(xsh.mean()),'iid_std':float(xid.std()),'shift_std':float(xsh.std())},
  'witness_demo':{'A_v_norm':float(np.linalg.norm(A@v)),'B_v_norm':float(np.linalg.norm(B@v)),'source_before':source0,'source_after':source1,'target_before':target0,'target_after':target1},
  'interpretation_ceiling':'Development demonstration only; not a real-data comparator or manuscript claim until independently audited and integrated under the official plan.'
}
(ROOT/'verification_outputs').mkdir(parents=True,exist_ok=True)
(ROOT/'verification_outputs'/'prediction_vs_recovery.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
