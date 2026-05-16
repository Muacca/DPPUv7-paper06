"""
KK Photon Field Strength
========================

Constructs the modified field strength F̃_{jk} and the O(A¹) connection
perturbation ω⁽¹⁾ for Kaluza-Klein photon reduction.

Physical background
-------------------
5D coframe:
    e^j = r₀ σ^j  (j = 0,1,2 : spatial S³/Nil³/T³ directions)
    e^3 = L (dτ + A_k σ^k)  (S¹ fiber with KK gauge field A_k)

Exterior derivative of e^3:
    de^3 = L [dA_k ∧ σ^k + A_k dσ^k]

The A_k dσ^k term generates topology-dependent corrections:
    T³  (dσ^k = 0)                      : F̃ = F
    Nil³ (dσ² = σ⁰∧σ¹)                  : F̃₀₁ = F₀₁ + A₂
    S³  (dσ^k ≠ 0)                      : F̃_{jk} = F_{jk} + 2ε_{jkl}A_l
    Sol³(dσ¹ = σ⁰∧σ¹, dσ² = -σ⁰∧σ²)     : F̃₀₁ = F₀₁ + A₁, F̃₀₂ = F₀₂ - A₂

The Koszul formula then gives the O(A¹) connection perturbation ω⁽¹⁾:
    ω⁽¹⁾(3, j, k) = +(L/2r₀²) F̃(j, k)
    ω⁽¹⁾(j, 3, k) = -(L/2r₀²) F̃(j, k)   where F̃ uses indices (j, k)
    ω⁽¹⁾(j, k, 3) = -(L/2r₀²) F̃(j, k)   where F̃ uses indices (j, k)
    ω⁽¹⁾(otherwise) = 0

This is the same structure for all topologies; only F̃ differs.

Author: Muacca
"""

from typing import Callable, Dict, List, Tuple

from sympy import Expr, Integer, expand
from sympy.tensor.array import MutableDenseNDimArray

from ..utils.epsilon import epsilon_3d


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------

_N_SPATIAL = 3   # number of spatial frame directions (0, 1, 2)
_FIBER_IDX = 3   # frame index of the S¹ KK fiber direction


# ---------------------------------------------------------------------------
# F_plain : standard gauge field strength
# ---------------------------------------------------------------------------

def make_F_plain(
    j: int,
    k: int,
    dA: List[List],
) -> Expr:
    """
    Standard (non-modified) field strength F_{jk} = dA[j][k] - dA[k][j].

    Args:
        j, k : spatial indices in {0, 1, 2}
        dA   : 3×3 array of derivative symbols dA[j][k] = ∂_j A_k

    Returns:
        SymPy expression for F_{jk}

    Notes:
        Antisymmetric by construction: F_{jk} = -F_{kj}.
    """
    if j == k:
        return Integer(0)
    return dA[j][k] - dA[k][j]


# ---------------------------------------------------------------------------
# topology-specific correction dicts
# ---------------------------------------------------------------------------

def s3_corrections(A: List) -> Dict[Tuple[int, int], Expr]:
    """
    Correction dict for S³ topology.

    F̃_{jk} = F_{jk} + 2ε_{jkl}A_l

    Keys are (j, k) with j < k (lower-triangle only).
    The function :func:`make_F_tilde` antisymmetrizes automatically.

    Args:
        A : list of 3 SymPy symbols [A₀, A₁, A₂]

    Returns:
        dict {(0,1): +2A₂,  (0,2): -2A₁,  (1,2): +2A₀}
    """
    return {
        (0, 1):  2 * A[2],
        (0, 2): -2 * A[1],
        (1, 2):  2 * A[0],
    }


def nil3_corrections(A: List) -> Dict[Tuple[int, int], Expr]:
    """
    Correction dict for Nil³ topology.

    F̃₀₁ = F₀₁ + A₂  (only the (0,1) direction is modified)
    All other F̃_{jk} = F_{jk}.

    Args:
        A : list of 3 SymPy symbols [A₀, A₁, A₂]

    Returns:
        dict {(0,1): +A₂}
    """
    return {(0, 1): A[2]}


def t3_corrections() -> Dict:
    """
    Correction dict for T³ topology.

    T³ is flat: dσ^k = 0, so no structural corrections.

    Returns:
        empty dict
    """
    return {}


def sol3_corrections(
    A: List,
    scale=1,
) -> Dict[Tuple[int, int], Expr]:
    """
    Correction dict for Sol³ topology.

    The left-invariant coframe satisfies
        d sigma^1 = + sigma^0 ^ sigma^1
        d sigma^2 = - sigma^0 ^ sigma^2

    so the KK field strength picks up the self-referential corrections
        F̃_01 = F_01 + scale * A_1
        F̃_02 = F_02 - scale * A_2

    The optional ``scale`` lets callers represent an inhomogeneous
    prefactor such as R_inv = 1 / R(z) without duplicating the sign/index
    pattern in multiple scripts.

    Args:
        A     : list of 3 SymPy symbols [A0, A1, A2]
        scale : multiplicative prefactor (default: 1)

    Returns:
        dict {(0,1): +scale*A1, (0,2): -scale*A2}
    """
    return {
        (0, 1):  scale * A[1],
        (0, 2): -scale * A[2],
    }


def sol3_inhomogeneous_corrections(
    A: List,
    R_inv,
) -> Dict[Tuple[int, int], Expr]:
    """
    Convenience wrapper for Sol³ corrections with an arbitrary 1/R profile.

    Args:
        A     : list of 3 SymPy symbols [A0, A1, A2]
        R_inv : symbolic prefactor representing 1 / R(z)

    Returns:
        dict {(0,1): +R_inv*A1, (0,2): -R_inv*A2}
    """
    return sol3_corrections(A, scale=R_inv)


# ---------------------------------------------------------------------------
# F̃ : modified field strength
# ---------------------------------------------------------------------------

def make_F_tilde(
    j: int,
    k: int,
    A: List,
    dA: List[List],
    corrections: Dict[Tuple[int, int], Expr],
) -> Expr:
    """
    Modified (covariant) field strength F̃_{jk} for a given topology.

    F̃_{jk} = F_{jk} + correction_{jk}(A)

    where correction_{jk} is topology-dependent and encodes the
    Maurer-Cartan equation of the spatial manifold (S³, Nil³, T³).

    Args:
        j, k        : spatial indices in {0, 1, 2}
        A           : list of 3 SymPy symbols
        dA          : 3×3 array of derivative symbols dA[j][k] = ∂_j A_k
        corrections : dict {(j,k): expr} for j < k only; lower-triangle.
                      Provide as :func:`s3_corrections`,
                      :func:`nil3_corrections`, or :func:`t3_corrections`.

    Returns:
        SymPy expression for F̃_{jk} (antisymmetric in j,k)
    """
    base = make_F_plain(j, k, dA)

    # Look up correction for (min, max) and antisymmetrize
    if j < k:
        corr = corrections.get((j, k), Integer(0))
    elif j > k:
        corr = -corrections.get((k, j), Integer(0))
    else:
        corr = Integer(0)

    return base + corr


# ---------------------------------------------------------------------------
# ω⁽¹⁾ : O(A¹) connection perturbation
# ---------------------------------------------------------------------------

def make_omega1(
    F_tilde_fn: Callable[[int, int], Expr],
    r0,
    L,
) -> Callable[[int, int, int], Expr]:
    """
    Build the O(A¹) connection perturbation ω⁽¹⁾ as a callable.

    The connection perturbation comes from the Koszul formula applied to
    the A-linear part of the structure constants C⁽¹⁾:
        C⁽¹⁾(3, j, k) = (L/r₀²) F̃(j, k)

    Koszul: ω⁽¹⁾_{abc} = (1/2)(C⁽¹⁾_{abc} + C⁽¹⁾_{cba} - C⁽¹⁾_{bac})

    Non-zero results:
        ω⁽¹⁾(3, j, k) = +(L/2r₀²) F̃(j, k)   j, k ∈ {0,1,2}
        ω⁽¹⁾(j, 3, k) = -(L/2r₀²) F̃(j, k)   j, k ∈ {0,1,2}
        ω⁽¹⁾(j, k, 3) = -(L/2r₀²) F̃(j, k)   j, k ∈ {0,1,2}
        ω⁽¹⁾(3, 3, k) = 0   [Koszul exact cancellation]
        ω⁽¹⁾(otherwise) = 0

    Args:
        F_tilde_fn : callable(j, k) → expr.  Typically a closure over A, dA:
                     ``lambda j, k: make_F_tilde(j, k, A, dA, corrections)``
        r0         : S³ radius (SymPy symbol or number)
        L          : S¹ length (SymPy symbol or number)

    Returns:
        Callable omega1(a, b, c) → SymPy expression.
    """
    prefactor = L / (2 * r0 ** 2)

    def omega1(a: int, b: int, c: int) -> Expr:
        if a == _FIBER_IDX and b < _N_SPATIAL and c < _N_SPATIAL:
            return prefactor * F_tilde_fn(b, c)
        elif a < _N_SPATIAL and b == _FIBER_IDX and c < _N_SPATIAL:
            return -prefactor * F_tilde_fn(a, c)
        elif a < _N_SPATIAL and b < _N_SPATIAL and c == _FIBER_IDX:
            return -prefactor * F_tilde_fn(a, b)
        else:
            # Includes (3,3,k), (k,3,3), (3,k,3): all zero by Koszul
            return Integer(0)

    return omega1


def omega1_to_array(
    omega1_fn: Callable[[int, int, int], Expr],
    dim: int = 4,
) -> MutableDenseNDimArray:
    """
    Convert an omega1 callable to a MutableDenseNDimArray.

    Used when passing ω⁽¹⁾ into the full-Riemann pipeline that
    requires array inputs (e.g. :func:`dppu.curvature.riemann.compute_riemann_tensor`).

    Args:
        omega1_fn : callable(a, b, c) from :func:`make_omega1`
        dim       : manifold dimension (default: 4 for 5D Kaluza-Klein)

    Returns:
        MutableDenseNDimArray of shape (dim, dim, dim)
    """
    arr = MutableDenseNDimArray.zeros(dim, dim, dim)
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                arr[a, b, c] = omega1_fn(a, b, c)
    return arr
