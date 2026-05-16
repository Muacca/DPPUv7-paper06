"""
Γ×Γ Shortcut for O(A²) Ricci Scalar
=====================================

Computes the O(A²) part of the 5D Ricci scalar using only the
O(A¹) connection perturbation ω⁽¹⁾, without computing the full
Riemann tensor.

Physical justification
----------------------
Expand the 5D connection as Ω = Ω⁽⁰⁾ + ω⁽¹⁾ + O(A²).
The Riemann tensor expands as:

    R^a_{bcd}[Ω] = R^(0)_{bcd} + R^(1)_{bcd} + R^(2)_{bcd} + O(A³)

where:
    R^(0): background curvature (constant; independent of A)
    R^(1) = 0: vanishes because the background satisfies the EOM
    R^(2) = ω⁽¹⁾ × ω⁽¹⁾ terms (the Γ×Γ shortcut)

Palatini identity:
    R^(2) = Γ²(ω⁽¹⁾) + total derivatives

On a compact space, total derivatives integrate to zero, so only the
Γ×Γ terms survive in the action integral.

The O(A²) Ricci scalar is therefore:

    R^{Γ²} = Σ_{a,b,e} [ω⁽¹⁾(a,a,e) ω⁽¹⁾(e,b,b) - ω⁽¹⁾(a,b,e) ω⁽¹⁾(e,b,a)]

This function avoids computing all 4⁴ = 256 Riemann components and
is 10-100× faster than the full Riemann approach.

Author: Muacca
"""

from typing import Callable

from sympy import Expr, Integer, expand, cancel


def gamma_gamma_ricci(
    omega1_fn: Callable[[int, int, int], Expr],
    dim: int = 4,
    simplify_result: bool = False,
) -> Expr:
    """
    Compute the O(A²) Ricci scalar via the Γ×Γ shortcut.

    Formula (Euclidean orthonormal frame, metric = δ_{ab}):

        R^{Γ²} = Σ_{a,b} R^{(2),a}_{bab}

        R^{(2),a}_{bcd} = Σ_e [ω⁽¹⁾(a,e,c) ω⁽¹⁾(e,b,d) - ω⁽¹⁾(a,e,d) ω⁽¹⁾(e,b,c)]

    For the KK photon with ω⁽¹⁾ from :func:`dppu.kk.field_strength.make_omega1`,
    this reproduces:
        - Maxwell term  -(L²/2r₀⁴) Σ F̃²
        - Chern-Simons term (topologically protected)
        - Mass² term from structure constant corrections

    Args:
        omega1_fn       : callable(a, b, c) → SymPy expr.
                          From :func:`dppu.kk.field_strength.make_omega1`.
        dim             : manifold dimension (default: 4)
        simplify_result : if True, apply cancel() to the result (slower)

    Returns:
        SymPy expression for R^{Γ²} at O(A²).
    """
    R_GG = Integer(0)

    for b in range(dim):
        d = b  # Euclidean identity metric: only diagonal b=d contributes
        for a in range(dim):
            for e in range(dim):
                t1 = omega1_fn(a, a, e) * omega1_fn(e, b, d)
                t2 = omega1_fn(a, d, e) * omega1_fn(e, b, a)
                R_GG += (t1 - t2)

    R_GG = expand(R_GG)

    if simplify_result:
        R_GG = cancel(R_GG)

    return R_GG
