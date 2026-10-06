# PPCO / ATI MASTER BULK 02 — INDEPENDENT AUDIT

Date: 2026-10-05

Final verdict: **PASS WITH OPERATOR-SCOPE NARROWING**

Return archive SHA-256:
`b31eef39a6a68ad59ae5977a42ba060a39581cd6a81846319f9b2ca7fa0bb6a3`

## Integrity
- ZIP CRC: PASS
- Returned files: 12
- Manifest-controlled payloads: 11
- Manifest hash mismatches: 0
- ATI module: `0.4.0rc5`
- Wheel SHA-256: `c072e707af6bdfbc556dfb5df8f25f86bba5a2fa4b0601f5370ec36831ef8c5b`
- Python: 3.10.6
- NumPy: 1.26.4
- SciPy: 1.14.1

## G09 production scaling
G09 independently reproduces as `PASS_MEASURED_RANGE`.

- 14 synthetic Gram cases.
- 20 frozen HCP directed operator pairs.
- All 20 real pairs: `NONEXACT_WITHIN_NUMERICAL_MODEL`.
- All 20 numerical missing dimensions equal the frozen exact-rank authority.
- Runner-reported full 20-pair elapsed time: 1.117996 s.
- Coarse enclosing-sampler peak: 65,122,304 bytes; this is NOT the maximum observed memory.
- Conservative maximum from the per-pair hash-controlled table: 70,316,032 bytes (67.06 MiB).
- Maximum observed per-pair working-set increase: 6,762,496 bytes (6.45 MiB).
- Maximum observed synthetic-grid working set: 81,879,040 bytes (78.09 MiB).
- No dense 902,629 x 902,629 projector claim is permitted.

The evidence supports measured sparse/Gram execution at the tested row dimensions and hardware only.

## G10 operator/preprocessing perturbation
G10 independently reproduces as `PASS_WITH_OPERATOR_SENSITIVITY_NARROW_TO_FROZEN_INSTANCE`.

- 18 frozen scenarios.
- 360 nonself pair evaluations.
- 90 perturbed-self controls.
- All 360 nonself evaluations remained nonexact.
- Zero nonself missing-dimension changes.
- Zero numerical UNKNOWN states.
- 85/85 self controls under non-rowspace-invariant perturbations became nonexact.
- 5/5 positive row-rescaling controls remained exact.

This is strong finite-instance evidence that exactness depends on the effective operator, not an atlas name alone.

## Scientific interpretation
The frozen hashed digital operator claims remain valid. The stress tests do not calibrate a distribution of biological registration error and do not prove that arbitrary real-world preprocessing obeys the same operator. Manuscript wording must therefore distinguish:

1. exact/no-go statements for frozen declared digital operators; from
2. conditional/approximate statements about real preprocessing pipelines.

No ABIDE/HCP outcome rescoring, PPCO fitting/retuning, CoRR access, subject replacement or parcel rescue occurred.

## Gate disposition
- G09: CLOSED PASS.
- G10: CLOSED PASS WITH REQUIRED CLAIM NARROWING.
- G11: DO_NOT_ADD retained.
- G12: DO_NOT_ADD retained.
- MASTER BULK 02: CLOSED PASS.
- Next: MASTER BULK 03 — release/provenance.

## Nonblocking erratum
The enclosing 20-pair memory sampler missed short transients captured by nested per-pair samplers. Public reporting uses the larger per-case maxima. See `G09_MEMORY_REPORTING_ERRATUM.md`.
