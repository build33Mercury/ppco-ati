# Exact alternate verifier protocol v0.3

The floating-point ATI producer and its NumPy verifier are not enough to establish that the theorem fixtures themselves are correct.

For analytic/synthetic fixtures whose entries are integers or exact rationals, a separate SymPy path must:

1. reject binary floating-point inputs;
2. compute ranks exactly over the rationals;
3. if `row(B) subseteq row(A)`, construct an exact rational `M` and verify `B = M A` symbolically;
4. otherwise construct `v in ker(A)` and verify `A v = 0`, `B v != 0` symbolically;
5. never use this exact-rational mode to imply exactness for real atlas operators that arise from floating interpolation unless a defensible exact representation exists.

Release gate: every analytic fixture used to justify theorem behavior must pass both the numerical producer/verifier path and the exact rational path where the fixture is rational.
