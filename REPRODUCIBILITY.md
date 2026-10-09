# Reproducibility and verification

This repository provides a separately inspectable reference implementation for conditional linear-operator identifiability analysis. It does **not** contain the complete historical HCP/ABIDE preprocessing, model input, and participant-level benchmark replay package.

## Software checks

From an isolated copy of the source directory:

```bash
python -m pip install '.[test]'
python -m pytest -q tests
python examples/example_certificate.py
```

Optional synthetic demonstrations (run in a writable working copy because the scripts write generated JSON results into the ignored `verification_outputs/` directory):

```bash
python scripts/run_randomized_validation.py
python scripts/run_robustness_stress_validation.py
python scripts/run_prediction_vs_recovery_demo.py
python scripts/run_scaling_benchmark.py
```

Development, adversarial, synthetic and unit tests check properties of the **reference software**. Their results must not be represented as biological, cross-site or participant-level benchmark validation. The numerical checks and any example certificates are conditional on the supplied matrices and stated assumptions.

## Figures and accompanying data

Install the optional figure dependencies and render the plots from aggregate values:

```bash
python -m pip install '.[figures]'
python scripts/make_submission_figures.py
```

The program writes to `submission_figures/`; that generated folder is intentionally excluded from version control. Six unique aggregate source tables are under [`source_data/`](source_data/). They preserve unfavorable outcomes, but do not allow independent recreation of historical participant-level results.

## External-data limitations

The historical benchmarks require exact software, atlas derivations, masks, operator matrices, subject-level inputs and original acquisition/preprocessing provenance. These are **not packaged here** and have not been independently reproduced from raw ABIDE/HCP participant imaging through this software release. Missing authorities cannot be replaced by similarly named or newly derived datasets without establishing equivalence.

The finite five-atlas support conclusions use declared historical operator identities, summarized in `FROZEN_OPERATOR_IDENTITIES.csv`. That list is an identity ledger, **not** the operator bytes or independently verified subject-effective models. No CoRR participant outcome analysis is claimed.

The software may be cited independently of the manuscript; an eventual journal article requires its own data, ethics and availability disclosures.
