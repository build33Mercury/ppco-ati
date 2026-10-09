# Scientific scope and interpretation

ATI Reference assesses whether stated target quantities are identifiable from specified linear source observations. Certificates and uncertainty classifications are conditional on the exact digital operators, observation model, and declared numerical tolerances.

## Supported uses

- Under the stated fixed linear covariance model, universal exact recovery follows the **established** row-space containment criterion. The software checks finite operator instances; it does not claim to have discovered this theorem.
- The reference implementation distinguishes universal exact recovery, conditional or bounded scalar inference, non-identifiability, and numerical abstention. Model-assisted prediction is a separate question from recoverability.
- For the five historical HCP digital operator identities recorded here, the 20 directed inter-atlas transfers were evaluated as nonexact under unrestricted positive-semidefinite covariance, conditional on the finite reconstructed support and assumptions. The reported stacked support rank is 601; this is a fixed-instance result, not a biological or arbitrary-atlas universality claim.
- Prospective archival sufficiency is defined for a **declared target family** and requires preservation of the necessary linear span and full joint covariance, including cross-covariances.
- The Visium illustration concerns aggregation and loss of reverse identifiability; it does not establish broad transfer validity across molecular assays.

## Evidence boundaries

- Frozen HCP covariance results do not support a general PPCO advantage. The ABIDE comparison is sensitive to site weighting and composition, and later subject-mask analyses do not constitute untouched external validation.
- Geometry/performance associations not externally replicated must not be described as validated.
- A numerical certificate establishes a fact about a supplied model, not the validity of the atlas, its preprocessing, or a participant's subject-effective observation operator.
- Synthetic tests, timing measurements and aggregate figures do not replace independent raw-input processing of HCP/ABIDE cohorts or prove clinical utility.
- Provider-controlled atlas and participant data are not redistributed. The operator identity hashes are a provenance aid, not raw input or a rights license.

The accompanying manuscript must report primary, sensitivity and exploratory analyses distinctly, retain unfavorable results and disclose unavailable evidence. See [Reproducibility](REPRODUCIBILITY.md) and [Data availability](DATA_AVAILABILITY.md).
