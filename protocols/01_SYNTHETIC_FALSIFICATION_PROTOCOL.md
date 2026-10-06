# Synthetic theorem/certificate falsification protocol

Freeze before real-data registry use.

Test families: exactly nested and equal row spaces; nonnested hidden directions; rank deficiency; nearly nested and ill-conditioned operators; hard partitions; soft overlaps; signed operators; empty/constant rows; coordinate permutations; covariance/correlation counterexamples; operator perturbations; unit rescaling.

For each case freeze analytic ground truth where possible, producer output, independent verifier output, tolerance sweep, and expected PASS/FAIL/UNKNOWN_NUMERICAL. Include deliberately false certificates to prove the verifier rejects them.

Release gate: no real atlas certificate resource until the independent verifier accepts 100% of valid frozen fixtures and rejects 100% of deliberately corrupted fixtures.
