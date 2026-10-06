from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from .identified_sets import scalar_variance_identified_interval, identified_interval_record, verify_interval_endpoint


def main(argv=None):
    p=argparse.ArgumentParser(description='Compute an ATI scalar target-variance identified interval under K>=0, A K A^T=C_A, with optional trace(K)<=tau.')
    p.add_argument('source', help='source operator A .npy')
    p.add_argument('source_covariance', help='source covariance C_A .npy')
    p.add_argument('target_vector', help='scalar target row/vector b .npy')
    p.add_argument('--tau', type=float, default=None, help='optional trace(K) upper bound')
    p.add_argument('--rtol', type=float, default=1e-10)
    p.add_argument('--atol', type=float, default=1e-12)
    p.add_argument('--out', required=True)
    p.add_argument('--save-endpoints', default=None, help='optional .npz path for endpoint K matrices')
    args=p.parse_args(argv)
    A=np.load(args.source,allow_pickle=False)
    C=np.load(args.source_covariance,allow_pickle=False)
    b=np.load(args.target_vector,allow_pickle=False)
    r=scalar_variance_identified_interval(A,C,b,tau=args.tau,rtol=args.rtol,atol=args.atol,return_endpoint_witnesses=True)
    record=identified_interval_record(r)
    checks={}
    if 'lower_endpoint_K' in r:
        checks['lower_endpoint']=verify_interval_endpoint(A,C,b,r['lower_endpoint_K'],tau=args.tau)
    if 'upper_endpoint_K' in r:
        checks['upper_endpoint']=verify_interval_endpoint(A,C,b,r['upper_endpoint_K'],tau=args.tau)
    record['endpoint_verification']=checks
    Path(args.out).write_text(json.dumps(record,indent=2,sort_keys=True),encoding='utf-8')
    if args.save_endpoints:
        payload={}
        if 'lower_endpoint_K' in r: payload['lower_K']=r['lower_endpoint_K']
        if 'upper_endpoint_K' in r: payload['upper_K']=r['upper_endpoint_K']
        np.savez(args.save_endpoints,**payload)
    print(args.out)

if __name__=='__main__':
    main()
