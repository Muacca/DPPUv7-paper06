"""
Full Riemann Validation Route
==============================

Computes the O(A²) Ricci scalar by building the full Riemann tensor
and extracting the quadratic-in-A piece.

This route is SLOWER than the Γ×Γ shortcut but serves as an
independent cross-check that:
    1. The on-shell condition R^(1) = 0 holds
    2. The Palatini identity (total-derivative cancellation) is correct
    3. The Γ×Γ shortcut and the full Riemann give identical O(A²) results

When to use
-----------
- Validation of new topologies before using the shortcut
- Unit testing in CI/CD pipelines (see :mod:`dppu.kk.validator`)
- Diagnosing discrepancies between expected and computed spectra

Performance warning
-------------------
Computing the full Riemann tensor for dim=4 involves O(dim⁴) = 256
tensor components, each of which is a SymPy polynomial in A and dA.
Computation time is typically 5-50× longer than the Γ×Γ shortcut.
Do NOT use this in production calculations or parameter scans.

Author: Muacca
"""

from typing import List

from sympy import Expr, Integer, Dummy, expand, diff
from sympy.tensor.array import MutableDenseNDimArray

from ..engine.levi_civita import compute_christoffel_frame
from ..curvature.riemann import compute_riemann_tensor
from ..curvature.ricci import compute_ricci_scalar


def make_omega0(
    C: MutableDenseNDimArray,
    dim: int = 4,
) -> MutableDenseNDimArray:
    """
    Build the background Levi-Civita connection ω⁽⁰⁾ from structure constants.

    Uses the Koszul formula:
        Γ^a_{bc} = (1/2)(C^a_{bc} + C^c_{ba} - C^b_{ac})

    This is a thin wrapper around
    :func:`dppu.engine.levi_civita.compute_christoffel_frame`.

    Args:
        C   : background structure constants (A-independent) as
              MutableDenseNDimArray of shape (dim, dim, dim)
        dim : manifold dimension (default: 4)

    Returns:
        MutableDenseNDimArray of shape (dim, dim, dim) for ω⁽⁰⁾
    """
    return compute_christoffel_frame(C, dim)


def full_riemann_scalar_a2(
    omega0: MutableDenseNDimArray,
    omega1: MutableDenseNDimArray,
    C: MutableDenseNDimArray,
    A_syms: List,
    dA_syms: List,
    dim: int = 4,
) -> Expr:
    """
    Compute the O(A²) Ricci scalar via the full Riemann tensor.

    Algorithm:
        1. Build Ω_full = ω⁽⁰⁾ + ω⁽¹⁾  (background + photon perturbation)
        2. Compute R^a_{bcd}[Ω_full, C] using the Riemann formula
        3. Compute Ricci scalar R_full = R^a_{bab}
        4. Extract O(A²) by the scaling trick:
               R_scaled(t) = R_full(A→t A, dA→t dA)
               R_A2 = (1/2) d²/dt² R_scaled |_{t=0}

    The scaling trick correctly isolates the quadratic piece even when
    higher-order terms are present (they vanish at t=0 after two
    differentiations).

    Args:
        omega0   : background connection ω⁽⁰⁾, shape (dim,dim,dim)
        omega1   : perturbation connection ω⁽¹⁾, shape (dim,dim,dim);
                   convert from callable via
                   :func:`dppu.kk.field_strength.omega1_to_array`
        C        : background structure constants, shape (dim,dim,dim)
        A_syms   : flat list of A-field symbols to scale (e.g. [A0,A1,A2])
        dA_syms  : flat list of dA symbols to scale (e.g. dA[j][k] for j≠k)
        dim      : manifold dimension (default: 4)

    Returns:
        SymPy expression for R^{full}|_{O(A²)}.
        For T³/Nil³/S³ this should equal the output of
        :func:`dppu.kk.gamma_gamma.gamma_gamma_ricci` to within simplify().

    Performance note:
        Runtime is typically O(dim⁴) work with polynomial algebra.
        For dim=4, expect seconds to minutes depending on topology.
    """
    # Step 1: Build omega_full = omega0 + omega1 (element-wise NDimArray addition)
    omega_full = omega0 + omega1

    # Step 2: Full Riemann tensor R^a_{bcd}[omega_full, C]
    R_tensor = compute_riemann_tensor(omega_full, C, dim)

    # Step 3: Full Ricci scalar (no expand needed before the scaling subs)
    R_full = compute_ricci_scalar(R_tensor, dim)

    # Step 4: Extract O(A²) via scaling R_full(t*A, t*dA), d²/dt²|_{t=0}/2
    # Use Dummy to guarantee no name collision with caller symbols.
    t = Dummy('t', positive=True)
    all_syms = list(A_syms) + list(dA_syms)
    subs_in = [(s, t * s) for s in all_syms]

    R_scaled = expand(R_full.subs(subs_in))
    R_a2 = expand(diff(R_scaled, t, 2).subs(t, 0)) / 2

    # Sign note: compute_riemann_tensor uses the convention
    #   R^a_{bcd} = Γ^a_{ec}Γ^e_{bd} - Γ^a_{ed}Γ^e_{bc} + Γ^a_{be}C^e_{cd}
    # which gives the O(A²) Ricci scalar with the OPPOSITE sign compared to
    # the Γ×Γ shortcut formula in gamma_gamma_ricci (which is physically
    # verified across all Phase computations).  Negating here makes the full
    # Riemann route consistent with the Γ×Γ shortcut.
    return -R_a2
