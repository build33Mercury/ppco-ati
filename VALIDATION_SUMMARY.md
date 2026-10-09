# Software and numerical validation summary

- Unit tests: 51/51 completed successfully.
- Randomized numerical validation: 4,902/4,902 cases completed without recorded failure.
- Robustness stress testing: 11,101/11,101 cases completed without recorded failure.
- Exact support-rank verification: rank 601 for the specified five-atlas operator family.
- Measured scaling analysis: completed over the reported source/target row sizes and latent dimensions; conclusions are limited to the measured range.
- Operator perturbation analysis: 360/360 nonself comparisons remained non-identifiable with unchanged reference missing dimension; 85/85 noninvariant self perturbations became nonexact, while 5/5 positive row-rescaling controls remained exact.

The numerical implementation is intended to return an unresolved classification rather than force an exact/nonexact decision when the specified numerical criteria are not sufficiently stable.
