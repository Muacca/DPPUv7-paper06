"""
Levi-Civita Connection (Koszul Formula)
=======================================

Computes the Levi-Civita (torsion-free, metric-compatible) connection
in a frame basis using the general Koszul formula.

Mathematical Background:
    For a left-invariant frame on a Lie group with structure constants C^a_{bc},
    the Levi-Civita connection is given by the Koszul formula:

    Γ^a_{bc} = (1/2)(C^a_{bc} + C^c_{ba} - C^b_{ac})

    This formula works for:
    - Bi-invariant metrics (e.g., S³ = SU(2))
    - Non bi-invariant metrics (e.g., Nil³ = Heisenberg)

    The key requirement is that the frame is left-invariant, not that
    the metric is bi-invariant.

CONVENTIONS Compliance:
    - Structure constants: [E_b, E_c] = -C^a_{bc} E_a (CONVENTIONS 3.2)
    - Koszul formula: CONVENTIONS 5

Author: Muacca
"""

from typing import Any, Optional

from sympy import S, cancel
from sympy.tensor.array import MutableDenseNDimArray


def compute_christoffel_frame(
    C: MutableDenseNDimArray,
    dim: int,
    logger: Optional[Any] = None,
    metric: Optional[Any] = None,
    metric_inv: Optional[Any] = None,
) -> MutableDenseNDimArray:
    """
    Compute Levi-Civita connection coefficients in frame basis.

    Uses the general Koszul formula (CONVENTIONS 5):
        Γ^a_{bc} = (1/2)(C^a_{bc} + C^c_{ba} - C^b_{ac})

    This formula is valid for any left-invariant frame, including
    non-bi-invariant cases like Nil³ (Heisenberg).

    Args:
        C: Structure constants C^a_{bc} as MutableDenseNDimArray
        dim: Dimension (typically 4)
        logger: Optional logger for progress messages
        metric: Optional frame metric g_{ab} for index lowering
        metric_inv: Optional inverse frame metric g^{ab} for index raising

    Returns:
        Connection coefficients Γ^a_{bc} as MutableDenseNDimArray

    Notes:
        - The resulting connection is automatically torsion-free
        - For orthonormal frames, it is also metric-compatible
        - Non-zero components are simplified before storage
    """
    if logger:
        logger.info("Computing Frame Connection Γ^a_bc via general Koszul formula...")

    Gamma_LC = MutableDenseNDimArray.zeros(dim, dim, dim)

    def _lower_c(c_arr: MutableDenseNDimArray, a: int, b: int, c: int) -> Any:
        if metric is None:
            return c_arr[a, b, c]
        return sum(metric[a, d] * c_arr[d, b, c] for d in range(dim))

    # GENERAL Koszul formula:
    # Γ_{abc} = (1/2)(C_{cbe} - C_{bec} - C_{ecb})
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                val = S.Zero
                for e in range(dim):
                    gamma_low = (
                        _lower_c(C, c, b, e)
                        - _lower_c(C, b, e, c)
                        - _lower_c(C, e, c, b)
                    ) / S(2)

                    if metric_inv is not None:
                        val += metric_inv[a, e] * gamma_low
                    else:
                        if a == e:
                            val += gamma_low

                if val != S.Zero:
                    Gamma_LC[a, b, c] = cancel(val)

    if logger:
        # Log non-zero components
        nonzero_count = 0
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    if Gamma_LC[a, b, c] != S.Zero:
                        nonzero_count += 1
                        if nonzero_count <= 10:
                            logger.info(f"  Γ^{a}_{{{b}{c}}} = {Gamma_LC[a,b,c]}")
        if nonzero_count > 10:
            logger.info(f"  ... and {nonzero_count - 10} more non-zero components")
        logger.success("Frame Connection computed")

    return Gamma_LC


def koszul_formula_biinvariant(
    C: MutableDenseNDimArray,
    a: int, b: int, c: int
) -> Any:
    """
    Koszul formula for bi-invariant metrics (simplified).

    For bi-invariant metrics (e.g., SU(2)):
        Γ^a_{bc} = (1/2) C^a_{bc}

    This is a special case of the general formula when the metric
    is bi-invariant (ad-invariant).

    Args:
        C: Structure constants
        a, b, c: Indices

    Returns:
        Connection coefficient Γ^a_{bc}

    Notes:
        This is provided for reference; the general formula should
        be used in the main computation pipeline.
    """
    return C[a, b, c] / S(2)


def check_torsion_free(
    Gamma: MutableDenseNDimArray,
    C: MutableDenseNDimArray,
    dim: int,
    logger: Optional[Any] = None
) -> bool:
    """
    Verify that a connection is torsion-free.

    Torsion-free condition (frame basis):
        T^a_{bc} = Γ^a_{bc} - Γ^a_{cb} + C^a_{bc} = 0

    Args:
        Gamma: Connection coefficients Γ^a_{bc}
        C: Structure constants C^a_{bc}
        dim: Dimension
        logger: Optional logger

    Returns:
        True if connection is torsion-free
    """
    if logger:
        logger.info("Checking torsion-free condition...")

    violations = 0
    for a in range(dim):
        for b in range(dim):
            for c in range(b + 1, dim):  # Only check b < c
                torsion = Gamma[a, b, c] - Gamma[a, c, b] + C[a, b, c]
                torsion = cancel(torsion)
                if torsion != S.Zero:
                    violations += 1
                    if logger and violations <= 3:
                        logger.warning(f"  T^{a}_{{{b}{c}}} = {torsion}")

    if logger:
        if violations == 0:
            logger.success("Torsion-free condition: PASSED")
        else:
            logger.error(f"Torsion-free condition: FAILED ({violations} violations)")

    return violations == 0
