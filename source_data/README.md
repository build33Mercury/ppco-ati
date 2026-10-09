# Aggregate figure source data

The six files here are distinct aggregate source tables for the scientific results and figure reconstruction. They are **not participant-level data** or a complete benchmark reproduction package.

| File | Contents |
| --- | --- |
| `Benchmark_contrast_authority.csv` | Frozen ABIDE/HCP covariance benchmark summary |
| `Figure2_covariance_RFE.csv` | Method-level covariance relative Frobenius errors |
| `Figure3_overlap_minus_PPCO.csv` | Overlap-minus-PPCO contrasts with reported uncertainty |
| `Figure4_operator_perturbation_self_controls.csv` | Operator perturbation self-control outcomes |
| `Figure5_scaling_measurements.csv` | Recorded operator-analysis timing measurements |
| `Operator_perturbation_family_summary.csv` | Perturbation-family aggregate results |

`python scripts/make_submission_figures.py` reconstructs Figures 1–5 using these files. The plot script does not rerun historical participant-level analyses. Figure 1 is schematic and contains no participant data.

The negative HCP comparison and site-sensitive ABIDE contrast must not be omitted or recast as evidence of general PPCO superiority.
