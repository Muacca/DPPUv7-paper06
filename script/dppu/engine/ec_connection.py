"""
Einstein-Cartan Connection
==========================

Computes the Einstein-Cartan (EC) connection from Levi-Civita
connection and contortion.

Mathematical Background:
    The Einstein-Cartan connection is the unique metric-compatible
    connection with prescribed torsion:

    ω^a_{bc} = Γ^a_{bc} + K^a_{bc}

    Where:
    - Γ^a_{bc} is the Levi-Civita (torsion-free) connection
    - K^a_{bc} is the contortion derived from torsion

Properties:
    - Metric-compatible: ∇ω g = 0
    - Torsion: T^a_{bc} = ω^a_{bc} - ω^a_{cb} + C^a_{bc}

References:
    - Hehl et al. (1976): Rev. Mod. Phys. 48, 393

Author: Muacca
"""

from typing import Any, Optional

from sympy import S, cancel
from sympy.tensor.array import MutableDenseNDimArray


def compute_ec_connection(
    Gamma_LC: MutableDenseNDimArray,
    K: MutableDenseNDimArray,
    dim: int,
    logger: Optional[Any] = None
) -> MutableDenseNDimArray:
    """
    Compute Einstein-Cartan connection from LC connection and contortion.

    Formula:
        ω^a_{bc} = Γ^a_{bc} + K^a_{bc}

    Args:
        Gamma_LC: Levi-Civita connection Γ^a_{bc}
        K: Contortion tensor K^a_{bc}
        dim: Dimension (typically 4)
        logger: Optional logger for progress messages

    Returns:
        EC connection ω^a_{bc} as MutableDenseNDimArray

    Notes:
        - If K = 0, this returns a copy of Gamma_LC
        - The result is NOT torsion-free (unless K = 0)
    """
    if logger:
        logger.info("Computing EC Connection ω^a_bc = Γ^a_bc + K^a_bc...")

    Gamma_EC = MutableDenseNDimArray.zeros(dim, dim, dim)

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                val = Gamma_LC[a, b, c] + K[a, b, c]
                Gamma_EC[a, b, c] = cancel(val) if val != S.Zero else S.Zero

    if logger:
        # Count non-zero components
        nonzero_lc = sum(
            1 for a in range(dim)
            for b in range(dim)
            for c in range(dim)
            if Gamma_LC[a, b, c] != S.Zero
        )
        nonzero_ec = sum(
            1 for a in range(dim)
            for b in range(dim)
            for c in range(dim)
            if Gamma_EC[a, b, c] != S.Zero
        )
        logger.info(f"  LC components: {nonzero_lc}, EC components: {nonzero_ec}")
        logger.success("EC Connection computed")

    return Gamma_EC


def verify_ec_torsion(
    Gamma_EC: MutableDenseNDimArray,
    C: MutableDenseNDimArray,
    T_expected: MutableDenseNDimArray,
    dim: int,
    logger: Optional[Any] = None
) -> bool:
    """
    Verify that EC connection has the expected torsion.

    Torsion formula (frame basis):
        T^a_{bc} = ω^a_{bc} - ω^a_{cb} + C^a_{bc}

    Args:
        Gamma_EC: EC connection ω^a_{bc}
        C: Structure constants C^a_{bc}
        T_expected: Expected torsion tensor
        dim: Dimension
        logger: Optional logger

    Returns:
        True if torsion matches expected
    """
    from sympy import cancel

    if logger:
        logger.info("Verifying EC connection has correct torsion...")

    violations = 0
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                T_computed = (
                    Gamma_EC[a, b, c] - Gamma_EC[a, c, b] + C[a, b, c]
                )
                diff = cancel(T_computed - T_expected[a, b, c])
                if diff != S.Zero:
                    violations += 1
                    if logger and violations <= 3:
                        logger.warning(
                            f"  T^{a}_{{{b}{c}}} mismatch: "
                            f"computed={T_computed}, expected={T_expected[a,b,c]}"
                        )

    if logger:
        if violations == 0:
            logger.success("EC torsion: VERIFIED")
        else:
            logger.error(f"EC torsion: MISMATCH ({violations} violations)")

    return violations == 0


def decompose_connection(
    Gamma_EC: MutableDenseNDimArray,
    Gamma_LC: MutableDenseNDimArray,
    dim: int
) -> MutableDenseNDimArray:
    """
    Extract contortion from EC connection by subtracting LC connection.

    Formula:
        K^a_{bc} = ω^a_{bc} - Γ^a_{bc}

    This is the inverse of compute_ec_connection.

    Args:
        Gamma_EC: EC connection ω^a_{bc}
        Gamma_LC: Levi-Civita connection Γ^a_{bc}
        dim: Dimension

    Returns:
        Contortion tensor K^a_{bc}
    """
    K = MutableDenseNDimArray.zeros(dim, dim, dim)

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                K[a, b, c] = Gamma_EC[a, b, c] - Gamma_LC[a, b, c]

    return K
