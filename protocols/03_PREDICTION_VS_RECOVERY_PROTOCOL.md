# Prospective prediction-versus-recovery stress test v0.2

**Freeze date:** 1 October 2026. **Status:** outcome-blind synthetic development protocol. No HCP E2, CoRR, or ABIDE LOCK outcomes may be used to alter it.

## Frozen observation model

- Latent dimension: 2.
- Source operator: `A = [[1,0]]`.
- Target operator: `B = [[0,1]]`.
- Hidden direction: `v = [0,1]`, so `Av=0` and `Bv=1`.
- Latent covariance samples are diagonal: `K_i = diag(x_i, y_i)`.
- Source observation is scalar `C_A=x_i`; target is scalar `C_B=y_i`.

## Frozen data-generating populations

Seed: `20261001`.

Training: `n=500`; `x ~ Uniform(0.5,2.0)`; `y = 1.0 + 0.8*x + eps`, `eps ~ Normal(0,0.05)`.

IID test: `n=500`; same generator as training.

Shift test: `n=500`; same source distribution but `y = 3.0 - 0.8*x + eps`, `eps ~ Normal(0,0.05)`.

All generated y values must remain positive; otherwise the run fails rather than clipping.

## Frozen predictor and metrics

Predictor: ordinary least-squares linear regression with intercept, implemented directly with NumPy. No hyperparameter tuning.

Report: MAE and R^2 for IID and shift test; fitted intercept/slope; source-distribution summary; exact equality of the source observation before/after adding a fixed `t vv^T` witness example with `t=2`.

## Interpretation ceiling

The purpose is to show that strong in-distribution prediction can coexist with universal non-identifiability and can be distribution-dependent. It is not evidence that a particular real-world predictor will fail, nor a benchmark against Krakencoder/ViTAE-HGOT.
