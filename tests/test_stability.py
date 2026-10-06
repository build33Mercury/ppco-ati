import numpy as np
from ati.stability import perturbation_stability_sweep


def test_perturbation_sweep_is_deterministic():
    A=np.array([[1.,0.,0.],[0.,1.,0.]])
    B=np.array([[1.,1.,0.]])
    x=perturbation_stability_sweep(A,B,epsilons=[0,1e-9],trials=5,seed=7)
    y=perturbation_stability_sweep(A,B,epsilons=[0,1e-9],trials=5,seed=7)
    assert x==y
    assert x['rows'][0]['status_counts']['EXACT_WITHIN_NUMERICAL_MODEL']==1


def test_nonnested_remains_nonnested_under_tiny_target_jitter():
    A=np.array([[1.,0.,0.]])
    B=np.array([[0.,1.,0.]])
    x=perturbation_stability_sweep(A,B,epsilons=[0,1e-12],trials=5,seed=8,perturb_source=False,perturb_target=True)
    assert x['rows'][0]['status_counts']['NONEXACT_WITHIN_NUMERICAL_MODEL']==1
    assert x['rows'][1]['status_counts']['NONEXACT_WITHIN_NUMERICAL_MODEL']==5
