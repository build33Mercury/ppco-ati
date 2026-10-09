# ATI Reference

Reference software for **Atlas Transfer Identifiability (ATI)**, accompanying the research manuscript *Atlas Transfer Is an Inverse Problem: Identifiability and Partial Inference for Functional Connectomes* (manuscript in preparation for *Neuroinformatics*; not published).

Version **0.4.0rc7** is a code-only prerelease. It implements numerical certificates for a **specified linear observation model and digital source/target operators**. Exact covariance recovery is assessed through the classical row-space inclusion criterion; the package also implements numerical non-identifiability checks, conditional scalar identified intervals, and an unresolved/abstention outcome when numerical assumptions are insufficient. These mathematical principles are established results, not new general theorems claimed by this software.

## Install

Requires Python 3.10 or newer. From the source root:

```bash
python -m pip install .
python -m ati.cli --help
python -m ati.interval_cli --help
python examples/example_certificate.py
```

For tests and figure-generation dependencies:

```bash
python -m pip install '.[test]'
python -m pytest -q tests
python -m pip install '.[figures]'
python scripts/make_submission_figures.py
```

The two command-line programs are also installed as `ati-certify` and `ati-identify-scalar`.

## Repository contents

| Directory | Purpose |
| --- | --- |
| `ati/` | Numerical identifiability analysis, verification, identified intervals, and CLI programs |
| `tests/`, `examples/` | Software regression tests and a synthetic example |
| `protocols/`, `schema/` | Mathematical specifications and machine-readable output schemas |
| `source_data/` | Six distinct aggregate tables used for figures and benchmark summaries |
| `scripts/` | Figure generation, synthetic checks, and scaling demonstrations |
| `evidence/` | Small aggregate numerical-sensitivity and scaling records, including adverse findings |

See [Scientific scope](SCIENTIFIC_SCOPE.md) for interpretation boundaries, [Reproducibility](REPRODUCIBILITY.md) for exactly what can and cannot be independently checked with this repository, [Data availability](DATA_AVAILABILITY.md) for data-access restrictions, and [Source data](source_data/README.md) for the figure tables.

## Scientific interpretation

- A certificate is valid only for the **supplied digital operators, declared observation model, and numerical tolerances**. It does not validate preprocessing choices, atlas provenance, or subject-specific effective masks.
- Prediction under a model or prior does not establish universal information recovery. The row-space condition and related linear-algebra results are classical.
- Frozen HCP covariance benchmark results do **not** support a general PPCO advantage; the ABIDE contrast is sensitive to site composition. The manuscript separately distinguishes prespecified results from later sensitivities.
- Synthetic software tests and aggregate figure tables do **not** reproduce participant-level HCP or ABIDE analyses, establish biological validity, or support clinical use.

## Citation and licensing

Use [`CITATION.cff`](CITATION.cff) when citing the software. Cite the manuscript separately only after it is available in a citable form. This repository includes no manuscript DOI.

The authored software is distributed under the [BSD 3-Clause License](LICENSE), subject to the [license scope](LICENSE_SCOPE.md). Third-party datasets and atlases remain governed by their own terms. This software release neither bundles nor licenses their underlying image data.

Current source-file hashes are listed in [`SHA256SUMS.txt`](SHA256SUMS.txt). GitHub visibility and release availability must be confirmed at the actual repository URL before describing this source as publicly available.
