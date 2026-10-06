# ATI theory specification v0.1 (prospective)

## Observation model
For linear scientific representation operators A and B acting on a common latent vector x, the covariance model is C_A = A K A^T and C_B = B K B^T for latent K >= 0. This statement is valid only when the preprocessing/statistic commutes with the operator. Atlas-specific nonlinear preprocessing, censoring, correlation normalization, Fisher transforms, or inconsistent time windows require a separate observation-model audit.

## T1 Exact universal covariance recoverability
Foundational statement, not presumed novel: C_B is a deterministic function of C_A for every PSD K if row(B) is contained in row(A). Then B = M A and C_B = M C_A M^T. Conversely, if row(B) is not contained in row(A), choose v in ker(A) with Bv != 0; adding t vv^T to any PSD K leaves C_A unchanged and changes C_B.

Novelty policy: treat the elementary factorization as standard/antecedent mathematics unless the literature audit establishes a narrower novel specialization.

## T2 Predictor impossibility corollary
If row(B) is not contained in row(A), no deterministic predictor f(C_A) can be uniformly correct for C_B over unrestricted PSD K: K and K+t vv^T produce the same input to f and different targets. This does not say a learned predictor cannot perform well on a restricted population distribution.

## T3 Minimal augmentation
Let S_A=row(A), S_B=row(B). The minimum number of additional independent linear measurements required for exact B recovery is dim(S_A+S_B)-dim(S_A)=dim(S_B)-dim(S_A intersect S_B), numerically represented by rank(B(I-P_A)). A valid implementation must return a basis U_add and independently verify row(B) subset row([A;U_add]).

## T4 Universal representation for multiple targets
For target operators B_1,...,B_m, the smallest exact source row space is span_i row(B_i), dimension rank([B_1;...;B_m]). This is standard subspace algebra; project value is in validated construction, cost-aware design extensions, and scientific archive guidance.

## T5 Conditional trace-bounded surrogate error
With P the orthogonal projector onto row(A), R=B(I-P), and tr(K)<=tau, the proposed conservative bound is
||B K B^T - B P K P B^T||_2 <= tau [2||BP||_2 ||R||_2 + ||R||_2^2].
This requires independent proof audit, equality/counterexample search, unit sensitivity analysis, and usefulness checks. It is not a universal no-assumption bound.

## T6 Correlation observation model
Correlation discards marginal scales. The 2D counterexample A=I, B=[[1,1],[1,-1]], K1=I, K2=diag(4,1) gives identical source correlations and target off-diagonal correlations 0 versus 0.6. A general characterization remains a research target; do not claim a broad theorem before proof.

## Per-estimand rule
ATI class must be attached to an observation model and estimand. A globally non-identifiable target covariance does not imply every target contrast is unbounded.
