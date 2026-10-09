# ATI Reference — 0.4.0rc7

**Software-only prerelease reference implementation; not a complete neuroimaging benchmark.**

ATI Reference evaluates linear-operator atlas transfer: exact identifiability under a declared observation model, non-identifiability witnesses, finite-model numerical uncertainty (UNKNOWN), and conditional scalar identified intervals. Row-space containment and the linear algebra are classical; this code is a reference implementation and auditable claim-calibration workflow, not a new universal-recovery theorem.

## Install and independently check the reference code

Requires Python 3.10+ and NumPy 1.24+. In an isolated environment, install the local source (and, optionally, testing dependencies):

```console
python -m pip install .
python -m pip install '.[test]'
python -m pytest -q tests
python examples/example_certificate.py
```

The scripts `scripts/run_randomized_validation.py`, `scripts/run_robustness_stress_validation.py`, and `scripts/run_prediction_vs_recovery_demo.py` are synthetic/development tests. Run them only against a writable detached copy because some scripts write JSON receipts. Use `python scripts/run_randomized_validation.py --help` only if supported; inspect the source before altering any default seeds.

CLI entry points after installation are `ati-certify` and `ati-identify-scalar`. The examples show their expected schema and status conventions. A numerical EXACT/NO-GO result applies to the *supplied hashed digital operators*; it is not a claim that the same status holds for an unverified subject-effective operator.

## What this archive **does not** reproduce

It does not include participant-level data from ABIDE, HCP or CoRR; third-party atlas images; the proprietary or rights-unclear derived operator matrices; or the complete historical analysis environment and 670 external benchmark authorities. The archived PPCO/ATI benchmarking outcomes are **not reproducible merely by installing ATI Reference**. Separate historical audit receipts provide internal reconciliation, not an independent raw-imaging replay. No clinical validity or general performance advantage is claimed. The frozen benchmark's unfavorable HCP and site-sensitive ABIDE contrasts remain relevant.

A public software release can be made independently from permitted data redistribution, but the accompanying manuscript must disclose the missing full benchmark reproduction route honestly. The repository URL is not evidence of public accessibility. Confirm repository visibility and release URL separately before claiming public availability.

## Rights and provenance

Software code: BSD-3-Clause (see LICENSE). Third-party imaging and atlas inputs keep their original permissions; **this repository's license does not cover them**. Data acquisition and validation plans are in the accompanying detached R11 gate bundle, not in the released software itself.

## Status

The 0.4.0rc7 source descends from reviewed rc6 reference code plus packaging, documentation, and smoke-test repairs. It is **not** a raw-imaging reproduction package. Historical `v0.4.0-rc6` receipts refer to that earlier version and must not be represented as newly executed rc7 benchmarks. The repository's public accessibility and release state must be checked at the final cited URL.


## Integrity manifests

`PUBLIC_SAFE_SHA256SUMS.txt` verifies the current source-release tree excluding the manifest itself, any build output and Git internals. Older rc6 provenance files, if retained, describe rc6 only.
