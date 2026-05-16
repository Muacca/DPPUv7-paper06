"""
Ricci Tensor and Scalar
=======================

Computes the Ricci tensor and Ricci scalar from the Riemann tensor.

Mathematical Background:
    Ricci tensor (contraction of Riemann):
        R_{bd} = R^a_{bad}

    Ricci scalar (trace of Ricci tensor):
        R = R^a_a = η^{ab} R_{ab}

    For orthonormal frame with identity metric:
        R = R_{aa} = sum over diagonal

Properties:
    - Ricci tensor is symmetric for Levi-Civita connection
    - For EC connection, Ricci tensor may have antisymmetric part
    - Ricci scalar determines the Einstein-Hilbert action

Author: Muacca
"""

from typing import Any, Optional

from sympy import S, cancel, Matrix
from sympy.tensor.array import MutableDenseNDimArray


def compute_ricci_tensor(
    Riemann: MutableDenseNDimArray,
    dim: int,
    logger: Optional[Any] = None
) -> Matrix:
    """
    Compute Ricci tensor from Riemann tensor.

    Formula:
        R_{bd} = R^a_{bad}

    Args:
        Riemann: Riemann tensor R^a_{bcd} as MutableDenseNDimArray
        dim: Dimension (typically 4)
        logger: Optional logger for progress messages

    Returns:
        Ricci tensor R_{bd} as SymPy Matrix

    Notes:
        - Contracts over first and third indices
        - Result is simplified component-wise
    """
    if logger:
        logger.info("Computing Ricci Tensor R_bd = R^a_bad...")

    Ricci = Matrix.zeros(dim, dim)

    for b in range(dim):
        for d in range(dim):
            val = S.Zero
            for a in range(dim):
                val += Riemann[a, b, a, d]
            Ricci[b, d] = cancel(val)

    if logger:
        # Check symmetry
        is_symmetric = all(
            Ricci[i, j] == Ricci[j, i]
            for i in range(dim)
            for j in range(i + 1, dim)
        )
        logger.info(f"  Ricci tensor symmetric: {is_symmetric}")
        logger.success("Ricci Tensor computed")

    return Ricci


def compute_ricci_scalar(
    Riemann: MutableDenseNDimArray,
    dim: int,
    logger: Optional[Any] = None,
    metric_inv: Optional[Any] = None
) -> Any:
    """
    Compute Ricci scalar from Riemann tensor.

    Formula:
        R = η^{bd} R^a_{bad}

    For orthonormal frame (Euclidean):
        R = sum_a sum_b R^a_{bab}

    Args:
        Riemann: Riemann tensor R^a_{bcd}
        dim: Dimension
        logger: Optional logger
        metric_inv: Optional inverse metric η^{ab} for proper contraction

    Returns:
        Ricci scalar R as SymPy expression
    """
    if logger:
        logger.info("Computing Ricci Scalar R = R^a_bab...")

    R_scalar = S.Zero
    for b in range(dim):
        for d in range(dim):
            ricci_bd = sum(Riemann[a, b, a, d] for a in range(dim))
            if metric_inv is not None:
                R_scalar += metric_inv[b, d] * ricci_bd
            else:
                if b == d:
                    R_scalar += ricci_bd

    R_scalar = cancel(R_scalar)

    if logger:
        logger.info(f"  Ricci Scalar R = {R_scalar}")

    return R_scalar


def compute_ricci_scalar_from_tensor(
    Ricci: Matrix,
    metric_inv: Matrix,
    dim: int
) -> Any:
    """
    Compute Ricci scalar from Ricci tensor.

    Formula:
        R = η^{ab} R_{ab}

    For identity metric:
        R = R_{00} + R_{11} + R_{22} + R_{33}

    Args:
        Ricci: Ricci tensor R_{ab}
        metric_inv: Inverse metric η^{ab}
        dim: Dimension

    Returns:
        Ricci scalar R
    """
    R_scalar = S.Zero
    for a in range(dim):
        for b in range(dim):
            R_scalar += metric_inv[a, b] * Ricci[a, b]

    return cancel(R_scalar)


def decompose_ricci_tensor(
    Ricci: Matrix,
    dim: int
) -> tuple:
    """
    Decompose Ricci tensor into symmetric and antisymmetric parts.

    R_{ab} = R_{(ab)} + R_{[ab]}

    Where:
        R_{(ab)} = (1/2)(R_{ab} + R_{ba})  (symmetric)
        R_{[ab]} = (1/2)(R_{ab} - R_{ba})  (antisymmetric)

    For Levi-Civita connection, R_{[ab]} = 0.
    For EC connection, R_{[ab]} may be non-zero.

    Args:
        Ricci: Ricci tensor R_{ab}
        dim: Dimension

    Returns:
        Tuple of (symmetric_part, antisymmetric_part) as Matrices
    """
    from sympy import Rational

    R_sym = Matrix.zeros(dim, dim)
    R_antisym = Matrix.zeros(dim, dim)

    for a in range(dim):
        for b in range(dim):
            R_sym[a, b] = cancel(Rational(1, 2) * (Ricci[a, b] + Ricci[b, a]))
            R_antisym[a, b] = cancel(Rational(1, 2) * (Ricci[a, b] - Ricci[b, a]))

    return R_sym, R_antisym
