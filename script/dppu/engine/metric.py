"""
Frame Metric Utilities
======================

Provides utilities for working with the frame (orthonormal) metric.

In our framework, we work in an orthonormal frame basis where:
- The frame metric η_{ab} = diag(1, 1, 1, 1) (Euclidean signature)
- Indices are raised/lowered with this metric
- Metric compatibility: ∇η = 0

Author: Muacca
"""

from typing import Any, List, Optional, Tuple

from sympy import Matrix, S, cancel
from sympy.tensor.array import MutableDenseNDimArray


def create_frame_metric(dim: int = 4, signature: str = "euclidean") -> Matrix:
    """
    Create an orthonormal frame metric.

    Args:
        dim: Dimension (default 4)
        signature: Metric signature
            - "euclidean": diag(+1, +1, +1, +1)
            - "lorentzian": diag(-1, +1, +1, +1)

    Returns:
        SymPy Matrix representing η_{ab}

    Examples:
        >>> metric = create_frame_metric(4, "euclidean")
        >>> metric[0, 0]
        1
    """
    if signature == "euclidean":
        return Matrix.eye(dim)
    elif signature == "lorentzian":
        diag_entries = [-1] + [1] * (dim - 1)
        return Matrix.diag(*diag_entries)
    else:
        raise ValueError(f"Unknown signature: {signature}")


def verify_metric_compatibility(
    Gamma: MutableDenseNDimArray,
    metric: Matrix,
    dim: int,
    logger: Optional[Any] = None
) -> Tuple[bool, List[Tuple]]:
    """
    Verify metric compatibility: Gamma_abc + Gamma_bac = 0.

    Gamma_abc = eta_ad Gamma^d_bc  (lower the first index with the metric).

    For the Levi-Civita connection the lowered coefficients must be
    antisymmetric in the first two indices.  This formulation is correct
    for both Euclidean (eta=I) and Lorentzian (eta=diag(-,+,+,+)) metrics.

    Args:
        Gamma: Connection coefficients Gamma^a_{bc}
        metric: Frame metric eta_{ab}  (used for index lowering)
        dim: Dimension
        logger: Optional logger for output

    Returns:
        Tuple of (passed: bool, violations: list)
        - passed: True if all checks pass
        - violations: List of (a, b, c, value) tuples for failed checks
    """
    if logger:
        logger.info("Verifying metric compatibility (Gamma_abc + Gamma_bac = 0)...")

    violations = []

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                # Lower first index: Gamma_abc = eta_{ad} Gamma^d_{bc}
                gamma_abc = sum(metric[a, d] * Gamma[d, b, c] for d in range(dim))
                gamma_bac = sum(metric[b, d] * Gamma[d, a, c] for d in range(dim))
                check = cancel(gamma_abc + gamma_bac)
                if check != S.Zero:
                    violations.append((a, b, c, check))

    passed = len(violations) == 0

    if logger:
        if passed:
            logger.success("Metric compatibility: PASSED")
        else:
            logger.error(f"Metric compatibility: FAILED ({len(violations)} violations)")
            for a, b, c, val in violations[:5]:
                logger.error(f"  Gamma_{{{a}{b}{c}}} + Gamma_{{{b}{a}{c}}} = {val}")

    return passed, violations


def raise_index(
    tensor: MutableDenseNDimArray,
    metric_inv: Matrix,
    index_pos: int,
    dim: int,
) -> MutableDenseNDimArray:
    """
    Raise a tensor index using the inverse metric.

    NOT IMPLEMENTED for metric-aware (Lorentzian) tensors.
    Use explicit contraction with metric_inv in calling code.

    Raises:
        NotImplementedError: always — guards against silent Lorentzian sign errors.
    """
    raise NotImplementedError(
        "raise_index is not implemented for metric-aware tensors. "
        "Use explicit contraction with metric_inv."
    )


def lower_index(
    tensor: MutableDenseNDimArray,
    metric: Matrix,
    index_pos: int,
    dim: int,
) -> MutableDenseNDimArray:
    """
    Lower a tensor index using the metric.

    NOT IMPLEMENTED for metric-aware (Lorentzian) tensors.
    Use explicit contraction with metric in calling code.

    Raises:
        NotImplementedError: always — guards against silent Lorentzian sign errors.
    """
    raise NotImplementedError(
        "lower_index is not implemented for metric-aware tensors. "
        "Use explicit contraction with metric."
    )
