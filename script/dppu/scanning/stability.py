"""
Stability Analysis
==================

Analyzes the stability of effective potential V(r).

Stability Types:
- type-I: Local minimum with barrier (V increases from r=0)
- type-II: Local minimum, V decreases from r=0 (spontaneous nucleation)
- type-III: No local minimum in physical region
"""

import logging
import warnings
from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple
import numpy as np
from scipy.optimize import minimize, minimize_scalar

logger = logging.getLogger(__name__)


class StabilityType(Enum):
    """Stability classification."""
    TYPE_I = "type-I"    # Barrier at origin
    TYPE_II = "type-II"  # Rolling from origin
    TYPE_III = "type-III"  # Unstable


def analyze_stability(
    potential_func: Callable,
    V_param: float,
    eta: float,
    theta: float,
    r_min: float = 0.01,
    r_max: float = 1e6,
    boundary_threshold: float = 0.02
) -> Tuple[Optional[float], Optional[float], StabilityType]:
    """
    Analyze stability of potential configuration.

    Args:
        potential_func: V(r, V_param, eta, theta) callable
        V_param, eta, theta: Parameter values
        r_min, r_max: Search bounds
        boundary_threshold: Threshold for boundary detection

    Returns:
        (r0, delta_V, stability_type)
    """
    def v(r):
        return potential_func(r, V_param, eta, theta)

    try:
        res = minimize_scalar(v, bounds=(r_min, r_max), method='bounded')
    except (ValueError, ArithmeticError, FloatingPointError) as e:
        logger.warning("minimize_scalar failed: %s", e)
        return None, None, StabilityType.TYPE_III

    if not res.success:
        return None, None, StabilityType.TYPE_III

    r0 = res.x
    v_min = res.fun

    # Boundary check
    if r0 < boundary_threshold or r0 > r_max - boundary_threshold:
        return None, None, StabilityType.TYPE_III

    # Curvature check
    h = 1e-5
    try:
        d2v = (v(r0 + h) - 2 * v(r0) + v(r0 - h)) / h**2
    except (ValueError, ArithmeticError, FloatingPointError) as e:
        logger.warning("Curvature check failed at r0=%.4f: %s", r0, e)
        return None, None, StabilityType.TYPE_III

    if d2v <= 0:
        return None, None, StabilityType.TYPE_III

    # Determine barrier type
    r_test = r_min * 2
    dr = r_min * 0.1
    try:
        slope = (v(r_test + dr) - v(r_test)) / dr
    except (ValueError, ArithmeticError, FloatingPointError) as e:
        logger.warning("Slope evaluation failed at r=%.4f: %s", r_test, e)
        slope = 0

    r_samples = np.linspace(r_min, r0, 30)
    v_samples = [v(r) for r in r_samples]
    v_max_before = max(v_samples)

    if slope > 0:
        stability_type = StabilityType.TYPE_I
        delta_V = v_max_before - v_min
    else:
        stability_type = StabilityType.TYPE_II
        delta_V = v_samples[0] - v_min
        if delta_V < 0:
            delta_V = abs(v_min)

    return r0, delta_V, stability_type


def find_equilibrium_r(
    potential_func: Callable,
    r_lo: float = 0.5,
    r_hi: float = 20.0,
    n_grid: int = 50,
    boundary_tol: float = 0.05,
) -> Tuple[Optional[float], Optional[float]]:
    """Robustly find r₀ = argmin V(r) within [r_lo, r_hi].

    Strategy: coarse grid search → local refinement → boundary check.

    Args:
        potential_func: V(r) callable (single argument)
        r_lo, r_hi: Search bounds
        n_grid: Number of grid points for coarse search
        boundary_tol: Fraction of [r_lo, r_hi] used as boundary exclusion zone.
                      If r₀ lands within this fraction of either boundary, returns (None, None).

    Returns:
        (r0, v_min) on success, or (None, None) if r₀ is clamped to boundary.
    """
    r_grid = np.linspace(r_lo, r_hi, n_grid)
    try:
        v_grid = np.array([float(potential_func(r)) for r in r_grid])
    except Exception as e:
        logger.warning("find_equilibrium_r: grid evaluation failed: %s", e)
        return None, None

    if not np.all(np.isfinite(v_grid)):
        v_grid = np.where(np.isfinite(v_grid), v_grid, np.inf)

    idx_min = int(np.argmin(v_grid))
    r_approx = float(r_grid[idx_min])
    step = (r_hi - r_lo) / n_grid
    r_lo_local = max(r_lo, r_approx - 3 * step)
    r_hi_local = min(r_hi, r_approx + 3 * step)

    try:
        res = minimize_scalar(potential_func, bounds=(r_lo_local, r_hi_local), method='bounded')
        r0 = float(res.x)
        v_min = float(res.fun)
    except Exception as e:
        logger.warning("find_equilibrium_r: local refinement failed: %s", e)
        r0 = r_approx
        v_min = float(v_grid[idx_min])

    span = r_hi - r_lo
    if r0 < r_lo + boundary_tol * span or r0 > r_hi - boundary_tol * span:
        warnings.warn(
            f"find_equilibrium_r: r₀={r0:.4f} is near boundary [{r_lo},{r_hi}]. "
            "Returning (None, None).",
            stacklevel=2,
        )
        return None, None

    return r0, v_min


def scan_vacuum_3d(
    veff_func: Callable,
    alpha: float,
    theta: float = 0.0,
    r_starts: Optional[List[float]] = None,
    eta_starts: Optional[List[float]] = None,
    V_starts: Optional[List[float]] = None,
    r_min: float = 0.05,
    boundary_tol: float = 0.05,
    xatol: float = 1e-8,
    fatol: float = 1e-10,
    maxiter: int = 3000,
) -> Tuple[Optional[Dict], Optional[float]]:
    """
    Find the global minimum of V_eff_EC(r, eta, V, alpha, theta) by multi-start
    Nelder-Mead optimization in (r, eta, V) space.

    Suitable for EC vacua where the minimum lies away from the (eta=0, V=0) slice,
    unlike the 1D ``find_equilibrium_r``.

    Args:
        veff_func:  Callable f(r, eta, V, alpha, theta) → float
        alpha:      Weyl coupling value (typically < 0 for EC vacuum)
        theta:      Nieh-Yan coupling (default 0)
        r_starts:   Coarse r seeds (default: [0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0])
        eta_starts: Coarse eta seeds (default: [0.0, 0.5, 1.0, -0.5])
        V_starts:   Coarse V seeds   (default: [0.0, 0.3, 0.5])
        r_min:      Minimum r boundary (rejected if r0 < r_min)
        boundary_tol: Unused (kept for API consistency)
        xatol, fatol, maxiter: Nelder-Mead options

    Returns:
        (params_dict, v_min) on success, where params_dict has keys
        {'r', 'eta', 'V', 'alpha', 'theta'}.
        Returns (None, None) if no interior minimum found.
    """
    if r_starts is None:
        r_starts = [0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
    if eta_starts is None:
        eta_starts = [0.0, 0.5, 1.0, -0.5]
    if V_starts is None:
        V_starts = [0.0, 0.3, 0.5]

    def objective(x: np.ndarray) -> float:
        r0, e0, v0 = x
        if r0 <= r_min:
            return 1e50
        try:
            val = float(veff_func(r0, e0, v0, alpha, theta))
            return val if np.isfinite(val) else 1e50
        except Exception:
            return 1e50

    best_x: Optional[np.ndarray] = None
    best_val = np.inf

    for r0 in r_starts:
        for e0 in eta_starts:
            for v0 in V_starts:
                try:
                    res = minimize(
                        objective, [r0, e0, v0],
                        method='Nelder-Mead',
                        options={'xatol': xatol, 'fatol': fatol,
                                 'maxiter': maxiter},
                    )
                    if res.x[0] > r_min and res.fun < best_val:
                        best_val = float(res.fun)
                        best_x = res.x.copy()
                except Exception:
                    pass

    if best_x is None or best_val >= 1e40:
        return None, None

    params_dict: Dict = {
        'r':     float(best_x[0]),
        'eta':   float(best_x[1]),
        'V':     float(best_x[2]),
        'alpha': float(alpha),
        'theta': float(theta),
    }
    return params_dict, best_val
