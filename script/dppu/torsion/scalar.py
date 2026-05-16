"""
Torsion Scalar
==============

Computes the torsion scalar T = T_{abc} T^{abc}.

Physical Background:
    The torsion scalar appears in various gravity Lagrangians:
    - Teleparallel gravity: L ~ T
    - f(T) gravity: L = f(T)
    - Our EC + NY theory: contributes to Nieh-Yan term

For orthonormal frame with identity metric:
    T = sum_{a,b,c} T_{abc}^2

Author: Muacca
"""

from typing import Any, Optional

from sympy import S, cancel
from sympy.tensor.array import MutableDenseNDimArray


def compute_torsion_scalar(
    T: MutableDenseNDimArray,
    dim: int = 4,
    logger: Optional[Any] = None
) -> Any:
    """
    Compute torsion scalar T = T_{abc} T^{abc}.

    For orthonormal frame with identity metric η_{ab} = δ_{ab}:
        T = sum_{a,b,c} T_{abc} T_{abc}

    Args:
        T: Torsion tensor T_{abc} as MutableDenseNDimArray
        dim: Dimension (default 4)
        logger: Optional logger

    Returns:
        Torsion scalar T as SymPy expression
    """
    if logger:
        logger.info("Computing Torsion Scalar T = T_abc T^abc...")

    T_scalar = S.Zero
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                T_scalar += T[a, b, c] * T[a, b, c]

    T_scalar = cancel(T_scalar)

    if logger:
        logger.info(f"  Torsion Scalar T = {T_scalar}")

    return T_scalar


def compute_torsion_pseudoscalar(
    T: MutableDenseNDimArray,
    dim: int = 4
) -> Any:
    """
    Compute torsion pseudoscalar (Pontryagin-like term).

    Definition:
        T_pseudo = (1/4) ε^{abcd} T_{eab} T^e_{cd}

    This appears in the TT part of the Nieh-Yan term.

    Args:
        T: Torsion tensor T_{abc}
        dim: Dimension

    Returns:
        Torsion pseudoscalar as SymPy expression
    """
    from sympy import Rational
    from ..utils.epsilon import epsilon_4d

    T_pseudo = S.Zero
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    eps = epsilon_4d(a, b, c, d)
                    if eps != 0:
                        sum_term = S.Zero
                        for e in range(dim):
                            sum_term += T[e, a, b] * T[e, c, d]
                        T_pseudo += Rational(1, 4) * eps * sum_term

    return cancel(T_pseudo)


def decompose_torsion_scalar(
    T: MutableDenseNDimArray,
    dim: int = 4
) -> dict:
    """
    Decompose torsion scalar into contributions from different parts.

    For our M³×S¹ ansatz with T = T1 + T2:
        T = T1·T1 + 2·T1·T2 + T2·T2

    Where:
        T1·T1: Pure axial contribution
        T2·T2: Pure vector-trace contribution
        T1·T2: Cross term (often zero by symmetry)

    Args:
        T: Torsion tensor
        dim: Dimension

    Returns:
        Dictionary with 'axial', 'vector', 'cross' contributions
    """
    # This requires separating T into T1 and T2 first
    # For now, return total scalar only
    total = compute_torsion_scalar(T, dim)
    return {
        'total': total,
        # Detailed decomposition would require mode information
    }


def compute_squared_norms(
    T: MutableDenseNDimArray,
    dim: int = 4
) -> dict:
    """
    Compute various squared norms of the torsion tensor.

    Norms computed:
        1. ||T||^2 = T_{abc} T^{abc}  (total)
        2. ||T_antisym||^2 = T_{[abc]} T^{[abc]}  (axial part)
        3. ||T_trace||^2 = T_a T^a  (trace vector norm)

    Args:
        T: Torsion tensor
        dim: Dimension

    Returns:
        Dictionary with 'total', 'antisym', 'trace' norms
    """
    # Total norm
    total = compute_torsion_scalar(T, dim)

    # Trace vector norm
    trace = []
    for b in range(dim):
        val = S.Zero
        for a in range(dim):
            val += T[a, a, b]
        trace.append(val)

    trace_norm = S.Zero
    for b in range(dim):
        trace_norm += trace[b] * trace[b]
    trace_norm = cancel(trace_norm)

    return {
        'total': total,
        'trace_norm': trace_norm,
    }
