"""
Pontryagin Density
==================

P = <R, *R> = R_flat · (*R)_flat

Functions:
- compute_P_from_riemann: P from numeric R^{ab}_{cd} array
- get_riemann_numeric: Evaluate symbolic Riemann from engine data
"""

from typing import Dict

import numpy as np

from .hodge import compute_hodge_dual


def compute_P_from_riemann(R_arr: np.ndarray) -> float:
    """
    Compute Pontryagin density P = <R, *R>.

    P = R_flat · (*R)_flat  (flat inner product over all indices)

    Args:
        R_arr: numpy array (4,4,4,4) representing R^{ab}_{cd}

    Returns:
        Scalar P value
    """
    R_dual = compute_hodge_dual(R_arr)
    return float(np.dot(R_arr.ravel(), R_dual.ravel()))


def get_riemann_numeric(
    engine_data: Dict,
    r: float,
    eta: float,
    V: float,
    alpha: float = 0.0,
    theta: float = 0.0,
    L: float = 1.0,
    kappa: float = 1.0,
) -> np.ndarray:
    """
    Evaluate symbolic Riemann tensor from UnifiedEngine data at a parameter point.

    Args:
        engine_data: dict returned by UnifiedEngine (engine.data)
        r, eta, V, alpha, theta, L, kappa: Parameter values

    Returns:
        numpy array (4,4,4,4) of R^{ab}_{cd} evaluated numerically
    """
    params = engine_data['params']
    riemann_sym = engine_data['riemann_abcd']

    subs = {
        params['r']:        r,
        params['eta']:      eta,
        params['V']:        V,
        params['alpha']:    alpha,
        params['theta_NY']: theta,
        params['L']:        L,
        params['kappa']:    kappa,
    }

    dim = engine_data['dim']
    R_arr = np.zeros((dim, dim, dim, dim), dtype=float)

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    val = riemann_sym[a, b, c, d]
                    if val != 0:
                        R_arr[a, b, c, d] = float(val.subs(subs))

    return R_arr
