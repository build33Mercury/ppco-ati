# Scalar identified-set proof obligations v0.3

Before the scalar trace-bounded interval enters a manuscript theorem/proposition:

- independently prove uniqueness of the source-row-space block `K_SS` under `A K A^T = C_A` after row-space orthogonalization;
- prove infeasibility criteria for inconsistent `C_A`;
- prove the hidden-direction Cauchy-Schwarz bound and endpoint attainability;
- prove the unbounded upper endpoint without a trace bound whenever the scalar target sees `ker(A)`;
- verify all edge cases: rank-zero A, singular K_SS, a=0, h=0, redundant source rows, tau exactly minimal, tau below minimal, and numerical perturbation;
- search optimal-recovery / moment-problem / covariance-completion literature for equivalent formulas;
- classify final novelty honestly as standard, adaptation, or new specialization.

Numerical endpoint construction is evidence of implementation consistency, not a substitute for proof.
