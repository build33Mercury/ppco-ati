# Observation-model commutation checklist v0.1 — frozen before registry/data expansion

A certificate is valid only for a declared observation model. For each dataset or representation pair, answer every item before applying the covariance-congruence theorem.

1. **Latent sample identity** — Are A and B applied to the same underlying samples/time points/voxels/cells?
2. **Temporal/sample selection** — Is censoring/scrubbing identical before representation? Any atlas-specific dropped samples invalidates naive commutation.
3. **Nuisance regression** — Is regression performed in the common latent space with identical regressors? Parcel-specific regression after aggregation is a different operator/statistic.
4. **Filtering/detrending/centering** — Are linear temporal transforms identical and applied before A/B? Record ordering.
5. **Mask/domain** — Do A and B act on the same frozen latent domain? Record mask hashes and excluded locations.
6. **Registration/resampling** — Is the mapping to common coordinates fixed and representation-independent? Record interpolation, affine/warp identity, and uncertainty.
7. **Parcel/region weights** — Are rows of A/B exactly defined, including normalization, partial volume, overlap, holes, and empty rows?
8. **Statistic** — Covariance, correlation, precision, Fisher-z, t-statistic, mean/count, or another object must be named explicitly. Do not transfer a theorem between statistics without derivation.
9. **Nonlinear normalization** — Global-signal scaling, per-parcel variance normalization, library-size normalization, adaptive segmentation, thresholding, clipping, and nonlinear transforms require a separate derivation.
10. **Finite-sample estimator** — State whether the theorem concerns population moments or the estimator actually computed; include missingness and degrees-of-freedom conventions.
11. **Commutation test fixture** — On a frozen synthetic/raw sample, compute both paths and quantify discrepancy: statistic(Ax) versus the theorem-predicted operator action on the latent statistic.
12. **Classification** — EXACT_MODEL / APPROXIMATE_MODEL / INVALID_FOR_THIS_THEOREM. Approximate status must carry an error/tolerance justification.

Fail closed: INVALID_FOR_THIS_THEOREM representations may remain in descriptive/resource tables but cannot receive ATI exact/no-go covariance certificates under the invalid model.
