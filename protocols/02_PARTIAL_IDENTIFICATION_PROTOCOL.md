# Conditional identified-set protocol v0.3

Status: partially implemented for one scalar variance estimand; general SDP remains prospective.

## Implemented analytic primitive

For `K >= 0`, `A K A^T = C_A`, and scalar target variance `b K b^T`, the reference implementation first orthogonalizes the row space of `A` and checks that `C_A` is consistent with a PSD observed block.

If `b` lies completely in `row(A)`, the scalar variance is exactly identified.

If `b` has a component in `ker(A)` and no bounding constraint is declared, the upper endpoint is unbounded and the lower endpoint is 0; the implementation returns `UNBOUNDED_UPPER` rather than a fabricated finite interval.

Under a prospectively justified trace constraint `tr(K) <= tau`, let `K_SS` be the identified source-row-space block, `a^2 = b_S^T K_SS b_S`, `h = ||b_H||`, and `T = tau - tr(K_SS)`. For feasible `T >= 0`, the exact scalar interval is

`[(max(a - sqrt(T) h, 0))^2, (a + sqrt(T) h)^2]`.

The implementation constructs PSD endpoint matrices and independently checks source equality, PSD, trace, and target value. This is a narrow theorem candidate requiring a separate prior-art/proof audit before manuscript novelty claims.

## General identified sets remain prospective

For matrix-valued targets or arbitrary linear functionals, investigate

`F(C_A,tau) = {K >= 0: A K A^T = C_A, tr(K) <= tau}`

and optimize linear target functionals with a validated SDP solver. Return solver/version, primal/dual status, residuals, objective, tau source, and sensitivity.

Noisy `C_A`: equality to an estimated point is not automatically justified. Predefine a source-covariance uncertainty set and propagate it. Distinguish inconsistent input, numerical failure, and scientific non-identifiability.

Do not advertise convex spectral-norm maximization as a routine convex SDP.
