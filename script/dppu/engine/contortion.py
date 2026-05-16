"""
Contortion Tensor
=================

Computes the contortion tensor from the torsion tensor.

Mathematical Background:
    The contortion tensor K_{abc} relates the Einstein-Cartan connection
    to the Levi-Civita connection:

    ω^a_{bc} = Γ^a_{bc} + K^a_{bc}

    Where the contortion is derived from torsion:

    K_{abc} = (1/2)(T_{abc} + T_{bca} - T_{cab})

    This formula ensures that the EC connection has the correct torsion.

Properties:
    - K_{abc} = -K_{bac} (antisymmetric in first two indices)
    - The torsion of the EC connection equals the input torsion

Author: Muacca
"""

from typing import Any, Optional

from sympy import Rational, cancel, S
from sympy.tensor.array import MutableDenseNDimArray


def compute_contortion(
    T: MutableDenseNDimArray,
    dim: int,
    metric=None,
    logger: Optional[Any] = None,
) -> MutableDenseNDimArray:
    """
    Compute contortion tensor K^a_{bc} from torsion tensor T_{abc}.

    Formula (all-lowered intermediate):
        K_{abc} = (1/2)(T_{abc} + T_{bca} - T_{cab})

    If metric is provided (Lorentzian-aware path):
        K^a_{bc} = g^{ad} K_{dbc}   (raise first index with inverse metric)

    If metric is None (legacy Euclidean path):
        Returns K_{abc} directly (identity raise, correct for identity metric).

    Args:
        T:      Torsion tensor T_{abc} as MutableDenseNDimArray
        dim:    Dimension (typically 4)
        metric: Frame metric g_{ab} (SymPy Matrix); None for legacy Euclidean path
        logger: Optional logger for progress messages

    Returns:
        Contortion tensor K^a_{bc} as MutableDenseNDimArray

    Notes:
        - For zero torsion, returns zero contortion
        - Non-zero components are simplified before storage
        - The Euclidean identity-metric path is preserved when metric=None
    """
    if logger:
        logger.info("Computing Contortion K^a_bc from Torsion...")

    K_low = MutableDenseNDimArray.zeros(dim, dim, dim)

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                val = (T[a, b, c] + T[b, c, a] - T[c, a, b]) * Rational(1, 2)
                if val != S.Zero:
                    K_low[a, b, c] = cancel(val)

    if metric is None:
        # Legacy Euclidean path: K^a_bc = K_abc (identity metric)
        K = K_low
    else:
        # Metric-aware path: raise first index K^a_bc = g^{ad} K_dbc
        metric_inv = metric.inv()
        K = MutableDenseNDimArray.zeros(dim, dim, dim)
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    val = sum(metric_inv[a, d] * K_low[d, b, c] for d in range(dim))
                    if val != S.Zero:
                        K[a, b, c] = cancel(val)

    if logger:
        nonzero_count = sum(
            1 for a in range(dim)
            for b in range(dim)
            for c in range(dim)
            if K[a, b, c] != S.Zero
        )
        logger.info(f"  Non-zero contortion components: {nonzero_count}")
        logger.success("Contortion computed")

    return K


def verify_contortion_antisymmetry(
    K: MutableDenseNDimArray,
    dim: int,
    logger: Optional[Any] = None
) -> bool:
    """
    Verify that contortion is antisymmetric in first two indices.

    Property:
        K_{abc} + K_{bac} = 0

    This follows from the definition and antisymmetry of torsion.

    Args:
        K: Contortion tensor K_{abc}
        dim: Dimension
        logger: Optional logger

    Returns:
        True if antisymmetry holds
    """
    if logger:
        logger.info("Verifying contortion antisymmetry K_abc = -K_bac...")

    violations = 0
    for a in range(dim):
        for b in range(a + 1, dim):
            for c in range(dim):
                check = cancel(K[a, b, c] + K[b, a, c])
                if check != S.Zero:
                    violations += 1
                    if logger and violations <= 3:
                        logger.warning(
                            f"  K_{{{a}{b}{c}}} + K_{{{b}{a}{c}}} = {check}"
                        )

    if logger:
        if violations == 0:
            logger.success("Contortion antisymmetry: PASSED")
        else:
            logger.error(f"Contortion antisymmetry: FAILED ({violations} violations)")

    return violations == 0


def contortion_from_axial_vector(
    S_mu: list,
    dim: int
) -> MutableDenseNDimArray:
    """
    Compute contortion from an axial torsion vector.

    For totally antisymmetric torsion (axial part):
        T_{abc} = (2/r) * ε_{abc} * S^τ  (in frame basis)

    The corresponding contortion is:
        K_{abc} = (1/r) * ε_{abc} * S^τ

    Args:
        S_mu: Axial vector components [S^0, S^1, S^2, S^3]
        dim: Dimension

    Returns:
        Contortion tensor for axial torsion

    Notes:
        This is a helper for the axial-only (AX) mode.
    """
    from ..utils.epsilon import epsilon_3d

    K = MutableDenseNDimArray.zeros(dim, dim, dim)

    # For M³×S¹, the axial part contributes to spatial indices only
    for a in range(3):
        for b in range(3):
            for c in range(3):
                eps = epsilon_3d(a, b, c)
                if eps != 0:
                    # K_{abc} = ε_{abc} * (axial component in τ direction)
                    K[a, b, c] = eps * S_mu[3]

    return K
