# Scaling memory-sampling reporting note

Classification: NONBLOCKING REPORTING / MEASUREMENT-SEMANTIC ERRATUM

The scaling summary reports an enclosing-block sampled working-set peak of
`65122304` bytes.

The same run's per-pair CSV contains higher inner-sampler observations. This is possible because both samplers poll asynchronously at 30 ms intervals; a nested per-pair sampler may observe a short transient missed by the enclosing sampler.

For conservative public reporting, use raw maxima from the preserved per-case tables:

- real 20-pair maximum observed working set: **70,316,032 bytes** (67.06 MiB);
- real 20-pair maximum observed within-case increase: **6,762,496 bytes** (6.45 MiB);
- synthetic-grid maximum observed working set: **81,879,040 bytes** (78.09 MiB);
- synthetic-grid maximum observed within-case increase: **14,856,192 bytes** (14.17 MiB).

The reported elapsed time and all exact/nonexact/missing-dimension results are unaffected.

No rerun is required because the raw per-case observations are present and hash-controlled. Future release/manuscript text must not quote the lower enclosing-sampler peak as the maximum observed memory.
