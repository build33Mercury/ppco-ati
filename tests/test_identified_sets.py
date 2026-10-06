import numpy as np
from ati.identified_sets import scalar_variance_identified_interval, verify_interval_endpoint


def test_exact_scalar_variance_when_target_in_source_rowspace():
    A=np.eye(2); C=np.array([[2.,.3],[.3,1.]])
    b=np.array([1.,2.])
    r=scalar_variance_identified_interval(A,C,b)
    truth=float(b@C@b)
    assert r['status']=='EXACT_WITHIN_NUMERICAL_MODEL'
    assert np.isclose(r['lower'],truth) and np.isclose(r['upper'],truth)


def test_unbounded_hidden_scalar_variance():
    A=np.array([[1.,0.]])
    C=np.array([[1.]])
    b=np.array([0.,1.])
    r=scalar_variance_identified_interval(A,C,b)
    assert r['status']=='UNBOUNDED_UPPER'
    assert r['lower']==0.0 and r['upper'] is None and r['upper_unbounded']
    assert verify_interval_endpoint(A,C,b,r['lower_endpoint_K'])['verified']
    assert np.isclose(verify_interval_endpoint(A,C,b,r['lower_endpoint_K'])['target_variance'],0.0,atol=1e-10)


def test_trace_bounded_interval_attains_both_endpoints():
    A=np.array([[1.,0.]])
    C=np.array([[1.]])
    b=np.array([1.,1.])
    r=scalar_variance_identified_interval(A,C,b,tau=2.0)
    assert r['status']=='BOUNDED_TRACE_CONSTRAINED'
    assert np.isclose(r['lower'],0.0,atol=1e-10)
    assert np.isclose(r['upper'],4.0,atol=1e-10)
    lo=verify_interval_endpoint(A,C,b,r['lower_endpoint_K'],tau=2.0)
    hi=verify_interval_endpoint(A,C,b,r['upper_endpoint_K'],tau=2.0)
    assert lo['verified'] and hi['verified']
    assert np.isclose(lo['target_variance'],r['lower'],atol=1e-9)
    assert np.isclose(hi['target_variance'],r['upper'],atol=1e-9)


def test_trace_bound_infeasible():
    A=np.array([[1.,0.]])
    C=np.array([[2.]])
    b=np.array([0.,1.])
    r=scalar_variance_identified_interval(A,C,b,tau=1.0)
    assert r['status']=='INFEASIBLE_TRACE_BOUND'


def test_source_covariance_infeasible_for_dependent_rows():
    A=np.array([[1.,0.],[2.,0.]])
    C=np.eye(2)
    b=np.array([1.,0.])
    r=scalar_variance_identified_interval(A,C,b,tau=10.)
    assert r['status']=='INFEASIBLE_SOURCE_COVARIANCE'


def test_random_feasible_completions_stay_inside_trace_interval():
    rng=np.random.default_rng(20261001)
    A=np.array([[1.,0.,0.],[0.,1.,0.]])
    C=np.array([[1.5,.2],[.2,.8]])
    b=np.array([.7,-.4,1.2])
    tau=3.0
    r=scalar_variance_identified_interval(A,C,b,tau=tau)
    assert r['status']=='BOUNDED_TRACE_CONSTRAINED'
    Ls=np.linalg.cholesky(C)
    slack=tau-np.trace(C)
    for _ in range(200):
        z=rng.normal(size=(1,2))
        zn=np.linalg.norm(z)
        if zn:
            z=z/zn*np.sqrt(slack)*rng.uniform(0,1)
        F=np.vstack([Ls,z])
        K=F@F.T
        val=float(b@K@b)
        assert r['lower']-1e-9 <= val <= r['upper']+1e-9
