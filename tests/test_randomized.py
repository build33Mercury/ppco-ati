import numpy as np
from ati import exact_recovery_analysis, minimal_augmentation, universal_representation, trace_bounded_projection_error_bound

RNG=np.random.default_rng(20261001)


def test_randomized_exact_nested_100():
    for _ in range(100):
        n=8; a=4; b=3
        A=RNG.normal(size=(a,n)); M=RNG.normal(size=(b,a)); B=M@A
        assert exact_recovery_analysis(A,B)['numerical_status']=='EXACT_WITHIN_NUMERICAL_MODEL'


def test_randomized_nonnested_100():
    for _ in range(100):
        n=8; a=3
        Q,_=np.linalg.qr(RNG.normal(size=(n,n)))
        A=Q[:,:a].T
        B=Q[:,a:a+2].T
        r=exact_recovery_analysis(A,B)
        assert r['numerical_status']=='NONEXACT_WITHIN_NUMERICAL_MODEL'
        assert r['witness'] is not None


def test_randomized_augmentation_100():
    for _ in range(100):
        n=8; A=RNG.normal(size=(3,n)); B=RNG.normal(size=(4,n))
        x=minimal_augmentation(A,B)
        assert x['verification_status']=='EXACT_WITHIN_NUMERICAL_MODEL'


def test_randomized_universal_50():
    for _ in range(50):
        ops=[RNG.normal(size=(2,7)), RNG.normal(size=(3,7)), RNG.normal(size=(1,7))]
        x=universal_representation(ops)
        assert all(s=='EXACT_WITHIN_NUMERICAL_MODEL' for s in x['target_verification_statuses'])


def test_randomized_trace_bound_250():
    for _ in range(250):
        n=7
        A=RNG.normal(size=(3,n)); B=RNG.normal(size=(4,n)); X=RNG.normal(size=(n,n)); K=X@X.T
        tau=float(np.trace(K))
        out=trace_bounded_projection_error_bound(A,B,tau)
        Q,_,vh=np.linalg.svd(A,full_matrices=False)
        rank=np.linalg.matrix_rank(A)
        basis=vh[:rank]
        P=basis.T@basis
        lhs=float(np.linalg.norm(B@K@B.T-B@P@K@P@B.T,2))
        assert lhs <= out['bound'] + 1e-8*max(1.0,out['bound'])
