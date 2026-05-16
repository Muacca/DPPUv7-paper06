"""
KK Effective Action Coefficient Extractor
==========================================

Extracts Maxwell, Chern-Simons, and mass² coefficients from the O(A²)
Ricci scalar returned by the Γ×Γ shortcut or the full Riemann route.

Structure of the KK photon effective action
-------------------------------------------
After integrating the 5D action over the internal manifold and keeping
only the O(A²) terms, the result has the form:

    S_eff ∝ ∫ d⁴x [ Σ_{i<j} k_ij (dA_{ij})² + k_CS ε^{ijk} A_i ∂_j A_k
                    + Σ_i k_i A_i² ]

These coefficients can be read off from the Ricci scalar:
    - Maxwell: coefficient of F²_{ij} = (dA_{ij})² per direction (i,j)
    - CS:      coefficient of F_{jk} × A_l (antisymmetric combination)
    - mass²:   coefficient of A_i² per spatial index i

API note
--------
``extract_maxwell`` and ``extract_mass`` return dicts so that anisotropic
spaces (e.g. Bianchi-IX squashing relevant for DPPUv3+) are handled
exactly.  An empty dict means the corresponding term is absent.  For
isotropic topologies (T³, S³) all dict values are equal; for Nil³ the
mass dict has only the massive direction(s).

Author: Muacca
"""

from typing import Dict, List, Tuple

from sympy import Expr, Integer, expand


def extract_maxwell(
    R_scalar: Expr,
    dA: List[List],
    n_spatial: int = 3,
) -> Dict[Tuple[int, int], Expr]:
    """
    Extract Maxwell term coefficients from the Ricci scalar.

    Scans all independent F²_{ij} components (i < j) and returns a dict
    mapping each (i, j) pair to its coefficient.  An empty dict means no
    Maxwell term is present.

    For isotropic spaces (T³, Nil³, S³ with uniform radius) all values
    are equal.  For anisotropic metrics (Bianchi-IX squashing) values
    differ between directions.

    Args:
        R_scalar  : SymPy expression (R^{Γ²} from gamma_gamma_ricci)
        dA        : n×n array of derivative symbols dA[j][k] = ∂_j A_k
        n_spatial : number of spatial directions (default: 3)

    Returns:
        Dict[(i, j) → Expr] for each (i < j) with a non-zero coefficient.
        Empty dict when no Maxwell term is found.
    """
    R_exp = expand(R_scalar)
    result: Dict[Tuple[int, int], Expr] = {}
    for i in range(n_spatial):
        for j in range(i + 1, n_spatial):
            coeff = R_exp.coeff(dA[i][j], 2)
            if coeff.is_zero is not True:
                result[(i, j)] = coeff
    return result


def extract_mass(
    R_scalar: Expr,
    A: List,
    dA: List[List],
) -> Dict[int, Expr]:
    """
    Extract photon mass² coefficients from the Ricci scalar.

    Returns a dict mapping each spatial index i to its A_i² coefficient.
    An empty dict means no mass term.

    Algorithm:
        1. Zero out all off-diagonal dA terms in one substitution pass.
        2. Extract A_i² coefficient for every spatial index i.
        3. Include only non-zero entries in the result dict.

    For isotropic spaces (T³, S³) all values are equal; for anisotropic
    spaces (Nil³: one massive axis, squashed S³) values differ.

    Args:
        R_scalar : SymPy expression
        A        : list of spatial SymPy symbols [A₀, A₁, A₂]
        dA       : n×n array of derivative symbols

    Returns:
        Dict[int → Expr] for each index with a non-zero coefficient.
        Empty dict when no mass term (e.g. T³).
    """
    n = len(A)
    zero_map = {
        dA[j][k]: Integer(0)
        for j in range(n)
        for k in range(n)
        if j != k
    }
    R_no_dA = expand(R_scalar).subs(zero_map)
    result: Dict[int, Expr] = {}
    for i in range(n):
        coeff = R_no_dA.coeff(A[i] ** 2)
        if coeff.is_zero is not True:
            result[i] = coeff
    return result


def extract_cs(
    R_scalar: Expr,
    A: List,
    dA: List[List],
) -> Expr | None:
    """
    Extract the Chern-Simons coefficient from the Ricci scalar.

    The CS term has the form k_CS × A₀ × (∂₁A₂ - ∂₂A₁) = k_CS × A₀ × F₁₂
    (and cyclic permutations, all with equal coefficient by S³ symmetry).

    Extraction reads the coefficient of A[0] × dA[1][2] directly,
    which equals k_CS (not k_CS × 2 because only one of the two antisymmetric
    terms is extracted).

    Args:
        R_scalar : SymPy expression
        A        : list of 3 SymPy symbols
        dA       : 3×3 array of derivative symbols

    Returns:
        SymPy expression for k_CS, or None if no CS term.
    """
    R_exp = expand(R_scalar)
    coeff = R_exp.coeff(A[0]).coeff(dA[1][2])
    return None if coeff.is_zero is True else coeff


def extract_all(
    R_scalar: Expr,
    A: List,
    dA: List[List],
) -> Dict[str, object]:
    """
    Extract all KK photon effective action coefficients.

    Convenience wrapper that calls all three extractors and returns
    a summary dict.

    Args:
        R_scalar : SymPy expression (R^{Γ²} or full Riemann O(A²))
        A        : list [A₀, A₁, A₂]
        dA       : 3×3 derivative symbol array

    Returns:
        dict with keys:
            'maxwell' : Dict[(i,j) → Expr] per direction.  Empty if none.
                        All values equal for isotropic spaces.
            'mass'    : Dict[int → Expr] per spatial index.  Empty if none.
                        e.g. {} for T³, {2: coeff} for Nil³,
                        {0: c, 1: c, 2: c} for S³.
            'cs'      : CS coefficient Expr (S³: -(2L²/r₀⁴)); None if absent.
    """
    return {
        'maxwell': extract_maxwell(R_scalar, dA),
        'mass':    extract_mass(R_scalar, A, dA),
        'cs':      extract_cs(R_scalar, A, dA),
    }
