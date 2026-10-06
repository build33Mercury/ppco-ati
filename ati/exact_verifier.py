from __future__ import annotations

from typing import Any, Dict, Iterable


def _sympy():
    try:
        import sympy as sp
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("sympy is required for exact rational verification") from exc
    return sp


def _rational_matrix(X: Iterable[Iterable[Any]]):
    sp = _sympy()
    rows = [list(r) for r in X]
    if not rows:
        return sp.Matrix([])
    width = len(rows[0])
    if any(len(r) != width for r in rows):
        raise ValueError("ragged matrix")
    def q(v):
        if isinstance(v, str):
            return sp.Rational(v)
        if isinstance(v, (int, sp.Integer, sp.Rational)):
            return sp.Rational(v)
        # Exact verification should never silently rationalize binary floats.
        if isinstance(v, float):
            raise TypeError("exact rational verifier rejects float inputs; use ints/strings/Rational")
        return sp.Rational(v)
    return sp.Matrix([[q(v) for v in r] for r in rows])


def exact_rational_classify(A, B) -> Dict[str, Any]:
    """
    Exact rational row-space containment test for analytic/synthetic fixtures.

    This is deliberately separate from the floating-point producer. It is not
    intended for resampled real atlas operators unless those operators have a
    defensible exact rational representation.
    """
    sp = _sympy()
    A = _rational_matrix(A)
    B = _rational_matrix(B)
    if A.cols != B.cols:
        raise ValueError("A and B must act on the same latent dimension")

    rank_a = int(A.rank())
    rank_stack = int(A.col_join(B).rank())
    exact = rank_stack == rank_a

    result: Dict[str, Any] = {
        "mode": "EXACT_RATIONAL",
        "rank_A": rank_a,
        "rank_B": int(B.rank()),
        "rank_stack": rank_stack,
        "rowspace_contained": bool(exact),
    }

    if exact:
        # Solve A.T * c = b.T independently for every target row.
        coeff_rows = []
        for i in range(B.rows):
            b = B.row(i).T
            solset = sp.linsolve((A.T, b))
            if solset is sp.EmptySet:
                raise AssertionError("rank test and exact solve disagreed")
            sol = next(iter(solset))
            # Free parameters can exist when A has dependent rows. Choose 0.
            free = sorted(set().union(*(x.free_symbols for x in sol)), key=lambda s: s.name)
            repl = {s: sp.Integer(0) for s in free}
            coeff_rows.append([sp.simplify(x.subs(repl)) for x in sol])
        M = sp.Matrix(coeff_rows)
        if M * A != B:
            raise AssertionError("exact factorization construction failed")
        result["exact_map_M"] = [[str(M[i, j]) for j in range(M.cols)] for i in range(M.rows)]
        result["verification"] = "B_EQUALS_MA_EXACTLY"
        return result

    # Find an exact null-space vector visible to B.
    witness = None
    for v in A.nullspace():
        bv = B * v
        if any(x != 0 for x in bv):
            # Normalize is unnecessary for exact algebra; keep primitive-ish vector.
            denoms = [sp.denom(x) for x in v]
            lcm = sp.ilcm(*[int(d) for d in denoms]) if denoms else 1
            vi = [sp.simplify(x * lcm) for x in v]
            nums = [abs(int(x)) for x in vi if x != 0 and x.is_integer]
            if nums:
                from math import gcd
                from functools import reduce
                g = reduce(gcd, nums)
                if g > 1:
                    vi = [sp.simplify(x / g) for x in vi]
            v = sp.Matrix(vi)
            witness = v
            break
    if witness is None:
        raise AssertionError("noncontainment found but no target-visible null witness constructed")
    if A * witness != sp.zeros(A.rows, 1):
        raise AssertionError("constructed witness is not in ker(A)")
    if B * witness == sp.zeros(B.rows, 1):
        raise AssertionError("constructed witness is not visible to B")

    result["witness_v"] = [str(x) for x in witness]
    result["Bv"] = [str(x) for x in (B * witness)]
    result["verification"] = "A_V_EQUALS_ZERO_AND_B_V_NONZERO_EXACTLY"
    return result
