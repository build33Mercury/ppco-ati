import copy
import numpy as np
from ati import (
    exact_recovery_analysis,
    find_kernel_witness,
    minimal_augmentation,
    universal_representation,
    principal_angles,
    generate_certificate,
    trace_bounded_projection_error_bound,
    build_interoperability_graph,
)
from ati.verifier import verify_certificate


def test_exact_identity_selection():
    A=np.eye(4); B=np.array([[1,0,0,0],[0,0,1,0]],float)
    r=exact_recovery_analysis(A,B)
    assert r['numerical_status']=='EXACT_WITHIN_NUMERICAL_MODEL'
    assert verify_certificate(generate_certificate(A,B),A,B)['verified']


def test_exact_linear_combination():
    A=np.array([[1,0,0],[0,1,0]],float); B=np.array([[2,-3,0]],float)
    r=exact_recovery_analysis(A,B)
    assert r['numerical_status']=='EXACT_WITHIN_NUMERICAL_MODEL'
    assert r['factorization_residual_relative'] < 1e-12


def test_nonexact_witness():
    A=np.array([[1,0,0]],float); B=np.array([[0,1,0]],float)
    r=exact_recovery_analysis(A,B)
    assert r['numerical_status']=='NONEXACT_WITHIN_NUMERICAL_MODEL'
    w=find_kernel_witness(A,B)
    assert w is not None and w['A_v_norm2'] < 1e-12 and w['B_v_norm2'] > 0.99
    assert verify_certificate(generate_certificate(A,B),A,B)['verified']


def test_minimal_augmentation_one_dimension():
    A=np.array([[1,0,0]],float); B=np.array([[1,1,0]],float)
    x=minimal_augmentation(A,B)
    assert x['minimal_dimension_numerical']==1
    assert x['verification_status']=='EXACT_WITHIN_NUMERICAL_MODEL'


def test_universal_representation():
    B1=np.array([[1,0,0]],float); B2=np.array([[0,1,0]],float); B3=np.array([[1,1,0]],float)
    u=universal_representation([B1,B2,B3])
    assert u['dimension_numerical']==2
    assert all(s=='EXACT_WITHIN_NUMERICAL_MODEL' for s in u['target_verification_statuses'])


def test_principal_angles_perpendicular():
    A=np.array([[1,0]],float); B=np.array([[0,1]],float)
    assert np.allclose(principal_angles(A,B),[np.pi/2],atol=1e-12)


def test_prediction_recovery_counterexample():
    A=np.array([[1,0]],float); B=np.array([[0,1]],float)
    K0=np.eye(2); v=np.array([0.,1.]); K1=K0+7*np.outer(v,v)
    assert np.allclose(A@K0@A.T,A@K1@A.T)
    assert not np.allclose(B@K0@B.T,B@K1@B.T)


def test_correlation_scale_counterexample():
    def corr(C):
        d=np.sqrt(np.diag(C)); return C/np.outer(d,d)
    A=np.eye(2); B=np.array([[1,1],[1,-1]],float)
    K1=np.eye(2); K2=np.diag([4.,1.])
    assert np.allclose(corr(A@K1@A.T),corr(A@K2@A.T))
    c1=corr(B@K1@B.T)[0,1]; c2=corr(B@K2@B.T)[0,1]
    assert abs(c1) < 1e-12 and np.isclose(c2,0.6)


def test_trace_bound_is_valid_for_example():
    A=np.array([[1.,0.,0.],[0.,1.,0.]])
    B=np.array([[1.,2.,3.],[-1.,0.,1.]])
    X=np.array([[1.,.2,0.],[.2,.5,.1],[0.,.1,.3]])
    K=X@X.T
    tau=float(np.trace(K))
    out=trace_bounded_projection_error_bound(A,B,tau)
    P=np.diag([1.,1.,0.])
    lhs=float(np.linalg.norm(B@K@B.T-B@P@K@P@B.T,2))
    assert lhs <= out['bound'] + 1e-10


def test_corrupted_exact_map_is_rejected():
    A=np.eye(3); B=np.array([[1.,0.,0.]])
    cert=generate_certificate(A,B)
    cert_bad=copy.deepcopy(cert)
    cert_bad['exact_map']['M'][0][0]=0.0
    check=verify_certificate(cert_bad,A,B)
    assert not check['verified'] and 'exact_map_residual_too_large' in check['issues']


def test_corrupted_hash_is_rejected():
    A=np.eye(3); B=np.array([[0.,1.,0.]])
    cert=generate_certificate(A,B)
    cert['source']['operator_sha256']='0'*64
    assert not verify_certificate(cert,A,B)['verified']


def test_corrupted_witness_is_rejected():
    A=np.array([[1.,0.,0.]])
    B=np.array([[0.,1.,0.]])
    cert=generate_certificate(A,B)
    cert['no_go_witness']['v']=[1.,0.,0.]
    check=verify_certificate(cert,A,B)
    assert not check['verified']
    assert 'witness_not_in_kernel' in check['issues'] or 'witness_not_visible_to_target' in check['issues']


def test_graph_direction_and_equivalence():
    ops={
        'full':np.eye(3),
        'x':np.array([[1.,0.,0.]]),
        'x_scaled':np.array([[2.,0.,0.]]),
        'y':np.array([[0.,1.,0.]])
    }
    g=build_interoperability_graph(ops)
    assert ['full','x'] in g['edges']
    assert ['x','full'] not in g['edges']
    assert ['x','x_scaled'] in g['edges'] and ['x_scaled','x'] in g['edges']


def test_unknown_numerical_band():
    A=np.array([[1.,0.]])
    B=np.array([[1.,1e-10]])
    r=exact_recovery_analysis(A,B,rtol=1e-10,atol=1e-12,stability_factor=10.0)
    assert r['numerical_status']=='UNKNOWN_NUMERICAL'
    cert=generate_certificate(A,B,rtol=1e-10,atol=1e-12)
    check=verify_certificate(cert,A,B)
    assert check['verified'] and check['mode']=='unknown_abstention'
