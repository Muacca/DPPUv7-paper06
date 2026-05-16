"""
Lagrangian Construction
=======================

Einstein-Cartan + Nieh-Yan Lagrangian: L = R/(2κ²) + θ_NY × N
"""

from typing import Any, Optional
from sympy import cancel, S


def compute_lagrangian(
    ricci_scalar: Any,
    nieh_yan_density: Any,
    kappa: Any,
    theta_NY: Any,
    weyl_scalar: Any = S.Zero,
    alpha: Any = S.Zero,
    logger: Optional[Any] = None
) -> Any:
    """
    Compute the EC + NY + Weyl Lagrangian density.

    L = R/(2κ²) + θ_NY × N + α × C²

    Args:
        ricci_scalar: Ricci scalar R
        nieh_yan_density: Nieh-Yan density N
        kappa: Gravitational coupling κ
        theta_NY: Nieh-Yan coupling θ_NY
        weyl_scalar: Weyl scalar C² (default 0)
        alpha: Weyl coupling α (default 0)

    Returns:
        Lagrangian density L
    """
    if logger:
        logger.info("Constructing Lagrangian L = R/(2κ²) + θ_NY×N + αC²...")

    L = ricci_scalar / (2 * kappa**2) + theta_NY * nieh_yan_density + alpha * weyl_scalar
    L = cancel(L)

    if logger:
        logger.info(f"  L = {L}")

    return L


def compute_action(
    lagrangian: Any,
    volume: Any,
    logger: Optional[Any] = None
) -> Any:
    """
    Compute total action S = ∫ L × Vol.

    Args:
        lagrangian: Lagrangian density L
        volume: Total volume element

    Returns:
        Action S
    """
    if logger:
        logger.info("Computing Action S = L × Vol...")

    action = cancel(lagrangian * volume)

    if logger:
        logger.info(f"  S = {action}")

    return action
