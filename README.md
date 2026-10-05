# ATI Reference 0.4.0rc6

Reference implementation accompanying **Atlas Transfer Is an Inverse Problem: Identifiability and Partial Inference for Functional Connectomes**.

ATI is a claim-calibration framework for declared scientific representation operators. The software distinguishes exact algebraic support, non-identifiability, numerically unresolved cases, and explicitly model/prior-assisted inference. It is designed to prevent representation conversion or predictive performance from being mislabeled as information recovery.

## What is established mathematics

The generic row-space/range-inclusion criterion, null-space witness construction, PSD covariance fibre, generic convex/SDP partial-identification machinery, and unrestricted span-augmentation dimension are treated as established foundations and are not claimed as new mathematics.

## Current software authority

`0.4.0rc6` is a release-only successor to the numerically validated `0.4.0rc5` authority. The numerical implementation is unchanged from rc5.

Validation:
- unit tests: 51/51 PASS;
- randomized development audit: 4,902/4,902 PASS;
- adversarial stress audit: 11,101/11,101 PASS;
- exact five-atlas support-rank adjudication: PASS;
- Windows Python 3.10 clean-room validation of the rc5 numerical authority: PASS;
- Master Bulk 02 measured scaling: PASS within the measured range;
- Master Bulk 02 operator-perturbation gate: PASS with operator-scope narrowing.

## Claim ceilings

- Exact/no-go claims apply to the declared hashed digital operators and stated observation model.
- A shared atlas name does not imply operator identity after registration, interpolation, support, weighting, or preprocessing changes.
- Prediction is not information recovery.
- Covariance and correlation are separate statistic spaces.
- PPCO is model-identified; there is no global PPCO-superiority claim.
- Archive-recoverability claims are target-family-specific.
- Scaling claims are limited to measured workloads.
- The repository does not redistribute atlas image bytes or derived atlas operator matrices.

## Bulk 02 operator perturbation result

Across 18 predeclared scenarios:
- 360/360 perturbed nonself transfers remained nonexact;
- 0/360 changed the frozen missing dimension;
- 0 numerical UNKNOWN states occurred;
- 85/85 non-rowspace-invariant perturbed self-controls became nonexact;
- 5/5 positive row-rescaling controls remained exact.

These are finite-instance digital-operator stress tests, not a calibrated distribution of biological registration error.

## Data and atlas policy

This repository contains no participant-level ABIDE/HCP/CoRR data, no atlas NIfTI images, and no derived project operator matrices. It publishes only safe code, protocols, aggregate evidence, and frozen identity metadata.

## License

ATI Reference is distributed under the **BSD-3-Clause** license. Atlas images and any third-party source data are not covered by this software license.
