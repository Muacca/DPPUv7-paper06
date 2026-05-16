"""
KK Route Validator
==================

Cross-checks the Γ×Γ shortcut against the full Riemann route for the
O(A²) Ricci scalar.

If both routes return identical expressions (up to SymPy simplification),
the following are simultaneously confirmed:
    - Palatini identity holds for this topology
    - on-shell condition R^(1) = 0 is satisfied by the background
    - Γ×Γ shortcut can be trusted for production calculations

Intended use
------------
Run once per new topology (or once after a significant code change) to
establish confidence, then use the fast Γ×Γ route for all subsequent work.

In CI/CD pipelines:
    from dppu.kk.validator import validate_kk_routes
    validate_kk_routes('s3', r0_val=3, L_val=1)  # raises on mismatch

Author: Muacca
"""

from typing import Optional

from sympy import symbols, Symbol, simplify, expand, Integer
from sympy.tensor.array import MutableDenseNDimArray

from .field_strength import (
    s3_corrections, nil3_corrections, t3_corrections,
    make_F_tilde, make_omega1, omega1_to_array,
)
from .gamma_gamma import gamma_gamma_ricci
from .full_riemann import make_omega0, full_riemann_scalar_a2


# ---------------------------------------------------------------------------
# Topology dispatch table (single source of truth)
# ---------------------------------------------------------------------------

_CORRECTIONS_DISPATCH = {
    's3':   lambda A: s3_corrections(A),
    'nil3': lambda A: nil3_corrections(A),
    't3':   lambda _A: t3_corrections(),
}

_KNOWN_TOPOLOGIES = list(_CORRECTIONS_DISPATCH.keys())


def _validate_topology(topology: str) -> str:
    """Normalise topology string and raise if unknown."""
    topo = topology.lower()
    if topo not in _CORRECTIONS_DISPATCH:
        raise ValueError(
            f"Unknown topology '{topology}'. Choose one of: {_KNOWN_TOPOLOGIES}"
        )
    return topo


def _build_topology_C(topology: str, r0, dim: int = 4) -> MutableDenseNDimArray:
    """
    Build the background (A-independent) structure constants for a given topology.

    Args:
        topology : 's3', 'nil3', or 't3' (case-insensitive)
        r0       : sphere radius
        dim      : manifold dimension

    Returns:
        MutableDenseNDimArray C[a,b,c] of shape (dim, dim, dim)

    Notes:
        Sign convention follows the dppu library (CONVENTIONS 3.2):
            [E_b, E_c] = -C^a_{bc} E_a

        S³:   C^i_{jk} = (4/r₀) ε_{ijk}
        Nil³: C^2_{01} = -1/r₀,  C^2_{10} = +1/r₀
        T³:   all zero
    """
    from ..utils.epsilon import epsilon_3d

    topo = _validate_topology(topology)
    C = MutableDenseNDimArray.zeros(dim, dim, dim)

    if topo == 's3':
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    eps = epsilon_3d(i, j, k)
                    if eps != 0:
                        C[i, j, k] = Integer(4) * eps / r0

    elif topo == 'nil3':
        C[2, 0, 1] = Integer(-1) / r0
        C[2, 1, 0] = Integer(1) / r0

    # t3: all zero (already initialised)

    return C


def validate_kk_routes(
    topology: str,
    r0=None,
    L=None,
    dim: int = 4,
    raise_on_mismatch: bool = True,
    verbose: bool = True,
) -> bool:
    """
    Validate that the Γ×Γ shortcut matches the full Riemann O(A²) result.

    Builds both routes symbolically for the given topology and checks
    that their difference simplifies to zero.

    Args:
        topology         : 's3', 'nil3', or 't3' (case-insensitive)
        r0               : radius symbol or positive number.
                           If None, uses ``symbols('r0', positive=True)``
        L                : length symbol or number.
                           If None, uses ``symbols('L', positive=True)``
        dim              : manifold dimension (default: 4)
        raise_on_mismatch: if True, raises AssertionError on failure
        verbose          : if True, prints progress and result

    Returns:
        True if both routes agree, False otherwise.

    Raises:
        AssertionError: if raise_on_mismatch=True and routes disagree.

    Example::

        from dppu.kk.validator import validate_kk_routes
        validate_kk_routes('t3')     # fastest, good sanity check
        validate_kk_routes('nil3')   # ~seconds
        validate_kk_routes('s3')     # ~tens of seconds
    """
    if r0 is None:
        r0 = symbols('r0', positive=True)
    if L is None:
        L = symbols('L', positive=True)

    topo = _validate_topology(topology)

    if verbose:
        print(f"[KK Validator] topology={topology}, r0={r0}, L={L}")

    # ── Build symbols ──
    A = [Symbol(f'A{k}') for k in range(3)]
    dA = [[Symbol(f'dA{j}{k}') if j != k else Integer(0)
           for k in range(3)] for j in range(3)]
    dA_flat = [dA[j][k] for j in range(3) for k in range(3) if j != k]

    # ── Topology-dependent corrections (via dispatch table) ──
    corrections = _CORRECTIONS_DISPATCH[topo](A)

    # ── F̃ and ω⁽¹⁾ (shared by both routes) ──
    F_tilde_fn = lambda j, k: make_F_tilde(j, k, A, dA, corrections)
    omega1_fn  = make_omega1(F_tilde_fn, r0, L)

    # ── Route 1: Γ×Γ shortcut ──
    if verbose:
        print("  [Route 1] Γ×Γ shortcut ...")
    R_shortcut = gamma_gamma_ricci(omega1_fn, dim=dim)

    # ── Route 2: full Riemann ──
    if verbose:
        print("  [Route 2] Full Riemann (may be slow) ...")

    C_arr = _build_topology_C(topo, r0, dim=dim)
    omega0_arr = make_omega0(C_arr, dim=dim)
    omega1_arr = omega1_to_array(omega1_fn, dim=dim)

    R_full = full_riemann_scalar_a2(
        omega0_arr, omega1_arr, C_arr,
        A_syms=A, dA_syms=dA_flat, dim=dim,
    )

    # ── Compare ──
    diff_expr = simplify(expand(R_shortcut - R_full))
    match = diff_expr.is_zero is True

    if verbose:
        status = "PASS ✓" if match else "FAIL ✗"
        print(f"  [Result] shortcut - full_riemann = {diff_expr}")
        print(f"  [Verdict] {status}")

    if not match and raise_on_mismatch:
        raise AssertionError(
            f"KK route mismatch for topology='{topology}':\n"
            f"  shortcut - full_riemann = {diff_expr}"
        )

    return match
