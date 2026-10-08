# Atlas Transfer Identifiability (ATI)

Reference software for **Atlas Transfer Is an Inverse Problem: Identifiability and Partial Inference for Functional Connectomes**, a research manuscript in preparation.

ATI evaluates claims about transfer between declared linear representations. Given specified source and target operators, it distinguishes universally exact recovery from non-identifiability and model-assisted prediction. Its certificates are conditional on the supplied operators and observation model. The general row-space criterion and associated inverse-problem mathematics are established results; this project applies them to specified digital atlas operators and supplies numerical verification tools.

## Install and inspect

Python 3.10 or newer is required.

```bash
python -m pip install -e ".[exact,test]"
python -m pytest -q
python examples/example_certificate.py
ati-certify --help
ati-identify-scalar --help
```

The example uses small synthetic operators. It does not download imaging data or reproduce the paper's participant-level benchmark.

## Repository contents

| Location | Contents |
| --- | --- |
| `ati/` | Certificate construction, verification, numerical checks and command-line interfaces. |
| `examples/`, `tests/` | A synthetic example and software tests. |
| `schema/` | Machine-readable certificate and registry schemas. |
| `protocols/` | Formal definitions, numerical checks and analysis specifications. |
| `source_data/` | Aggregate values used for manuscript figures. |
| `evidence/` | Archived software validation and finite-operator diagnostics. |

The historical evidence files retain their original identifiers so they can be matched to frozen reports. Those identifiers are not a recommended naming convention for new studies.

## Scope and limitations

- Exact recovery claims concern the declared, hashed digital operators. An atlas name alone does not identify an operator after registration, interpolation, masking or weighting.
- Predictive accuracy does not imply recovery of information lost by an observation operator.
- PPCO is a model-assisted estimator. The reported ABIDE and HCP covariance comparisons do not establish general PPCO superiority. A later HCP masked-operator sensitivity was conducted after earlier results were known and is reported separately in the manuscript.
- The five-atlas result is a finite-instance statement about those operators and the specified observation model. It is not a claim of biological validation or universal transfer across atlases.

See [scientific claims](ADMISSIBLE_CLAIMS.md) and the [atlas data policy](ATLAS_DATA_POLICY.md) for the precise boundaries.

## Data and reproducibility

This repository contains the ATI reference library, tests, protocols, aggregate figure values and selected validation evidence. It **does not contain the complete ABIDE/HCP benchmark workflow or a clean-room, raw-input replay package**. Participant imaging, source atlas images and derived project operator matrices are not included. Eligible users must obtain provider data through the relevant access routes and comply with their terms. The seven identified ABIDE atlas resource files can be independently acquired from the pinned upstream repository; that availability does not establish redistribution rights for the atlas binaries.

The immutable software tag `v0.4.0-rc6` identifies the archived reference implementation. The current default branch may receive documentation corrections without changing that tag. The manuscript and supplementary information describe the numerical results, adverse findings and remaining reproducibility limits; this repository should not be cited as a complete reproduction of the participant-level study.

## Citation and license

Use [CITATION.cff](CITATION.cff) for the software citation. The manuscript has no publication DOI at present. The ATI software is BSD-3-Clause licensed; this license does not cover provider imaging or third-party atlases. See [LICENSE](LICENSE) and [third-party notices](THIRD_PARTY_NOTICES.md).
