# ATI Reference

Python reference implementation for linear-operator identifiability and conditional scalar inference in functional-connectome atlas transfer. The software accompanies *Atlas Transfer Is an Inverse Problem: Identifiability and Partial Inference for Functional Connectomes* (manuscript in preparation).

## Installation

Python 3.10 or newer is required.

```bash
python -m pip install .
python -m ati.cli --help
python -m ati.interval_cli --help
python examples/example_certificate.py
```

Run the test suite:

```bash
python -m pip install '.[test]'
python -m pytest -q tests
```

Regenerate the aggregate-data figures:

```bash
python -m pip install '.[figures]'
python scripts/make_submission_figures.py
```

Outputs are written to `submission_figures/` (not committed). Additional numerical checks are available in `scripts/`.

## Files

- `ati/`: numerical analysis, certificate verification, identified intervals, and command-line tools
- `tests/` and `examples/`: regression tests and synthetic examples
- `schema/`: machine-readable output specifications
- `scripts/`: numerical validation, benchmarking demonstrations, and figure generation
- `source_data/`: six distinct aggregate tables for the accompanying figures
- `FROZEN_OPERATOR_IDENTITIES.csv`: recorded identities and dimensions of historical digital operators, **not** the operator matrices themselves

## Reproducibility and scope

Certificates concern **supplied digital operators and stated numerical/model assumptions**. The example tests do not validate participant-level data processing. This repository contains no HCP/ABIDE participant imaging or third-party atlas images, nor the full historical benchmark inputs and preprocessing workflow. Consequently it does **not** independently reproduce the manuscript's participant-level comparisons or establish clinical validity. The aggregate tables retain unfavorable as well as favorable results.

## Citation and license

For software citation see `CITATION.cff`. The software code is under the [BSD 3-Clause License](LICENSE); third-party data and atlas licensing are separate and no third-party source data or imaging files are redistributed here.
