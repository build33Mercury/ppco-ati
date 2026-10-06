import pytest
from ati.exact_verifier import exact_rational_classify


def test_exact_rational_containment_and_map():
    A=[[1,0,0],[0,1,0]]
    B=[[2,-3,0]]
    r=exact_rational_classify(A,B)
    assert r['rowspace_contained']
    assert r['verification']=='B_EQUALS_MA_EXACTLY'
    assert r['exact_map_M']==[['2','-3']]


def test_exact_rational_witness():
    A=[[1,0,0]]
    B=[[0,1,0]]
    r=exact_rational_classify(A,B)
    assert not r['rowspace_contained']
    assert r['verification']=='A_V_EQUALS_ZERO_AND_B_V_NONZERO_EXACTLY'
    assert r['witness_v'] is not None


def test_exact_rational_rejects_floats():
    with pytest.raises(TypeError):
        exact_rational_classify([[1.0,0.0]],[[1.0,0.0]])
