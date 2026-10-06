import copy
import json
from pathlib import Path
import numpy as np
import pytest
import jsonschema

from ati import exact_recovery_analysis, generate_certificate, minimal_augmentation, universal_representation, gram_containment_analysis
from ati.verifier import verify_certificate
from ati.identified_sets import scalar_variance_identified_interval, verify_interval_endpoint, identified_interval_record


def test_tau_infeasible_even_when_target_exact():
    A=np.array([[1.,0.]])
    C=np.array([[1.]])
    b=np.array([1.,0.])
    r=scalar_variance_identified_interval(A,C,b,tau=.5)
    assert r['status']=='INFEASIBLE_TRACE_BOUND'

@pytest.mark.parametrize('tau',[np.nan,np.inf,-1.])
def test_invalid_tau_rejected_before_exact_shortcut(tau):
    A=np.array([[1.,0.]])
    C=np.array([[1.]])
    b=np.array([1.,0.])
    with pytest.raises(ValueError): scalar_variance_identified_interval(A,C,b,tau=tau)


def test_asymmetric_source_covariance_fails_closed():
    A=np.eye(2); C=np.array([[1.,.8],[.1,1.]])
    with pytest.raises(ValueError): scalar_variance_identified_interval(A,C,np.array([1.,0.]))


def test_tiny_hidden_direction_is_unknown_not_exact():
    A=np.array([[1.,0.]])
    C=np.array([[1.]])
    b=np.array([1.,1e-11])
    r=scalar_variance_identified_interval(A,C,b)
    assert r['status']=='UNKNOWN_NUMERICAL'


def test_verifier_ignores_malicious_infinite_producer_rtol():
    A=np.array([[1.,0.]])
    B=np.array([[0.,1.]])
    cert=generate_certificate(A,B)
    # Replace a valid no-go payload with a false exact payload and malicious producer policy.
    cert.pop('no_go_witness',None)
    cert['exact_map']={'M':[[0.0]]}
    cert['claim']['ati_class']='ATI0_EXACT_ALGEBRAIC_NUMERICALLY_SUPPORTED'
    cert['numerical_policy']['rtol']=1e-10  # schema-valid; verifier still uses its own policy
    check=verify_certificate(cert,A,B)
    assert not check['verified']
    assert 'independent_rowspace_test_not_exact' in check['issues']


def test_verifier_rejects_unreasonable_local_tolerance():
    A=np.eye(2); B=A.copy(); cert=generate_certificate(A,B)
    with pytest.raises(ValueError): verify_certificate(cert,A,B,rtol=1.0)


def test_augmentation_archive_requires_cross_covariance():
    A=np.array([[1.,0.]])
    B=np.array([[1.,1.]])
    x=minimal_augmentation(A,B)
    req=x['archive_requirement']
    assert req['cross_covariance_required']
    assert not req['separate_C_AA_and_C_UU_sufficient']
    assert req['incremental_unique_float64_values']==2  # C_AU one value + C_UU one value


def test_cross_covariance_counterexample():
    # Same source variance and augmentation variance, opposite cross-covariances -> different target variance.
    A=np.array([[1.,0.]])
    U=np.array([[0.,1.]])
    B=np.array([[1.,1.]])
    Kp=np.array([[1.,.8],[.8,1.]])
    Km=np.array([[1.,-.8],[-.8,1.]])
    assert np.allclose(A@Kp@A.T,A@Km@A.T)
    assert np.allclose(U@Kp@U.T,U@Km@U.T)
    assert not np.allclose(B@Kp@B.T,B@Km@B.T)


def test_universal_archive_byte_formula():
    ops=[np.array([[1.,0.,0.]]),np.array([[0.,1.,0.]])]
    x=universal_representation(ops)
    assert x['dimension_numerical']==2
    assert x['archive_requirement']['packed_float64_bytes_per_subject']==24


def test_gram_ratio_is_sqrt_of_energy_ratio_and_agrees_with_direct():
    rng=np.random.default_rng(44)
    A=rng.normal(size=(3,8)); B=rng.normal(size=(4,8))
    Gaa=A@A.T; Gbb=B@B.T; Gba=B@A.T
    g=gram_containment_analysis(Gaa,Gbb,Gba,source_shape=A.shape,target_shape=B.shape)
    d=exact_recovery_analysis(A,B)
    assert g['numerical_status']==d['numerical_status']
    assert np.isclose(g['direct_residual_relative']**2,g['residual_energy_ratio'],rtol=1e-8,atol=1e-12)
    assert np.isclose(g['direct_residual_relative'],d['rowspace_residual_relative'],rtol=1e-7,atol=1e-10)


def test_interval_endpoint_rejects_asymmetric_K():
    A=np.eye(2); C=np.eye(2); b=np.array([1.,0.]); K=np.array([[1.,1.],[0.,1.]])
    with pytest.raises(ValueError): verify_interval_endpoint(A,C,b,K)


def test_schema_migrations_validate():
    root=Path(__file__).resolve().parents[1]
    A=np.eye(3); B=np.array([[1.,0.,0.]])
    cert=generate_certificate(A,B)
    schema=json.loads((root/'schema'/'ati_certificate_v0_3.schema.json').read_text())
    jsonschema.validate(cert,schema)
    C=np.array([[1.]])
    rec=identified_interval_record(scalar_variance_identified_interval(np.array([[1.,0.]]),C,np.array([1.,1.]),tau=2.,return_endpoint_witnesses=False))
    ischema=json.loads((root/'schema'/'ati_scalar_interval_v0_2.schema.json').read_text())
    jsonschema.validate(rec,ischema)


def test_gram_requires_original_operator_shapes():
    A=np.array([[1.,0.]])
    B=np.array([[1.,1.]])
    g=gram_containment_analysis(A@A.T,B@B.T,B@A.T)
    assert g['numerical_status']=='UNKNOWN_NUMERICAL'
    assert g['reason']=='original_operator_shapes_required_for_shared_rank_policy'


def test_gram_contract_uses_large_latent_shape_without_dense_projector():
    # Gram blocks do not encode the latent column count. Supply a large original shape and
    # verify that rank policy uses sqrt(eigenvalues) plus the original m x n dimensions.
    n=200000
    Gaa=np.diag([1.0,1e-14])   # implicit source singular values [1, 1e-7]
    Gbb=np.diag([1.0,0.0])
    Gba=np.array([[1.0,0.0],[0.0,0.0]])
    g=gram_containment_analysis(Gaa,Gbb,Gba,source_shape=(2,n),target_shape=(2,n))
    assert g['numerical_status']=='EXACT_WITHIN_NUMERICAL_MODEL'
    assert g['source_rank_numerical']==1
    assert g['target_rank_numerical']==1



def test_gram_identical_self_identity_short_circuits_cancellation():
    from ati.core import gram_containment_analysis
    # Deliberately ill-scaled positive diagonal Gram to force a realistic large-n
    # numerical contract while preserving exact self identity.
    d=np.geomspace(1e-6,1e-2,116)
    G=np.diag(d)
    out=gram_containment_analysis(G,G,G,source_shape=(116,902629),target_shape=(116,902629))
    assert out["numerical_status"]=="EXACT_WITHIN_NUMERICAL_MODEL"
    assert out["reason"]=="identical_gram_self_identity"
    assert out["missing_dimension_numerical"]==0
    assert out["direct_residual_relative"]==0.0


def test_gram_identical_self_identity_does_not_relax_nearby_nonidentity():
    from ati.core import gram_containment_analysis
    G=np.eye(4)
    Gba=G.copy(); Gba[0,0]=0.0
    out=gram_containment_analysis(G,G,Gba,source_shape=(4,1000),target_shape=(4,1000))
    assert out["numerical_status"]=="NONEXACT_WITHIN_NUMERICAL_MODEL"
