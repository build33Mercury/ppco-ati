from __future__ import annotations
import json, platform, sys, time, tracemalloc
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
import numpy as np
from ati import generate_certificate
from ati.verifier import verify_certificate

rng=np.random.default_rng(20261001)
sizes=[32,64,128,256]
reps=3
rows=[]
for n in sizes:
    a=max(4,n//4); b=max(2,n//8)
    A=rng.normal(size=(a,n)); M=rng.normal(size=(b,a)); B=M@A
    times=[]; peaks=[]
    for _ in range(reps):
        tracemalloc.start(); t0=time.perf_counter()
        cert=generate_certificate(A,B)
        check=verify_certificate(cert,A,B)
        dt=time.perf_counter()-t0
        _,peak=tracemalloc.get_traced_memory(); tracemalloc.stop()
        if not check['verified']: raise RuntimeError('exact certificate verification failed')
        times.append(dt); peaks.append(peak)
    rows.append({'latent_n':n,'source_rows':a,'target_rows':b,'case':'nested_exact',
                 'seconds_median':float(np.median(times)),'seconds_min':float(min(times)),'seconds_max':float(max(times)),
                 'python_tracemalloc_peak_bytes_median':int(np.median(peaks))})

    # strong nonnested target from orthogonal complement when possible
    Q,_=np.linalg.qr(rng.normal(size=(n,n)))
    A2=Q[:,:a].T; B2=Q[:,a:a+b].T
    times=[]; peaks=[]
    for _ in range(reps):
        tracemalloc.start(); t0=time.perf_counter()
        cert=generate_certificate(A2,B2)
        check=verify_certificate(cert,A2,B2)
        dt=time.perf_counter()-t0
        _,peak=tracemalloc.get_traced_memory(); tracemalloc.stop()
        if not check['verified']: raise RuntimeError('witness certificate verification failed')
        times.append(dt); peaks.append(peak)
    rows.append({'latent_n':n,'source_rows':a,'target_rows':b,'case':'nonnested_witness',
                 'seconds_median':float(np.median(times)),'seconds_min':float(min(times)),'seconds_max':float(max(times)),
                 'python_tracemalloc_peak_bytes_median':int(np.median(peaks))})

result={'artifact':'ATI Reference synthetic scaling check','development_only':True,
        'repetitions':reps,'rows':rows,'python':sys.version,'numpy':np.__version__,'platform':platform.platform(),
        'warning':'Synthetic CPU development timing only; not a claim of atlas-scale production scalability. tracemalloc excludes some native BLAS allocations.'}
Path('verification_outputs').mkdir(parents=True,exist_ok=True)
Path('verification_outputs/scaling_check.json').write_text(json.dumps(result,indent=2,sort_keys=True),encoding='utf-8')
print(json.dumps(rows,indent=2))
