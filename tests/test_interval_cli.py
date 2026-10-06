import json, subprocess, sys
from pathlib import Path
import numpy as np


def test_interval_cli_module(tmp_path):
    A=np.array([[1.,0.]])
    C=np.array([[1.]])
    b=np.array([1.,1.])
    np.save(tmp_path/'A.npy',A); np.save(tmp_path/'C.npy',C); np.save(tmp_path/'b.npy',b)
    out=tmp_path/'interval.json'
    subprocess.check_call([sys.executable,'-m','ati.interval_cli',str(tmp_path/'A.npy'),str(tmp_path/'C.npy'),str(tmp_path/'b.npy'),'--tau','2','--out',str(out)])
    x=json.loads(out.read_text())
    assert x['status']=='BOUNDED_TRACE_CONSTRAINED'
    assert x['endpoint_verification']['lower_endpoint']['verified']
    assert x['endpoint_verification']['upper_endpoint']['verified']
