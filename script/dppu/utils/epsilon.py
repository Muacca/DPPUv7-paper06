"""
Epsilon Symbols and Epsilon Tensors (Levi-Civita)
==================================================

Provides:
  - epsilon_symbol: combinatorial Levi-Civita symbol (signature-independent)
  - epsilon_4d: backward-compatible alias (symbol only, not tensor)
  - epsilon_tensor_down: lower-index epsilon tensor (metric-dependent)
  - epsilon_tensor_up: upper-index epsilon tensor (metric-dependent)
  - metric_signature_sign: sign(det(metric))
  - raise_epsilon_index_sign: product of diagonal metric_inv elements

Conventions (DPPU-LZ v1.1):
  epsilon_symbol(0,1,2,3) = +1  (signature-independent)

  Lorentzian orthonormal frame diag(-,+,+,+):
    epsilon_tensor_down(eta, 0,1,2,3) = +1
    epsilon_tensor_up(eta_inv, 0,1,2,3) = -1

  Euclidean orthonormal frame diag(+,+,+,+):
    epsilon_tensor_down(delta, 0,1,2,3) = +1
    epsilon_tensor_up(delta, 0,1,2,3) = +1
"""

from typing import Tuple


def epsilon_symbol(a: int, b: int, c: int, d: int) -> int:
    """
    Combinatorial Levi-Civita symbol. Signature-independent.
    Symbol only, not tensor.

    epsilon_symbol(0,1,2,3) = +1
    Returns -1, 0, or +1.
    """
    indices = [a, b, c, d]
    if len(set(indices)) != 4:
        return 0
    inversions = 0
    for i in range(4):
        for j in range(i + 1, 4):
            if indices[i] > indices[j]:
                inversions += 1
    return 1 if inversions % 2 == 0 else -1


def epsilon_4d(mu: int, nu: int, rho: int, sigma: int) -> int:
    """
    4D Levi-Civita symbol for indices 0, 1, 2, 3.
    Symbol only, not tensor. Signature-independent.

    Backward-compatible alias for epsilon_symbol.
    epsilon_4d(0,1,2,3) = +1.

    Do NOT use this as a Lorentzian epsilon tensor.
    For tensor operations use epsilon_tensor_down / epsilon_tensor_up.
    """
    return epsilon_symbol(mu, nu, rho, sigma)


def epsilon_3d(i: int, j: int, k: int) -> int:
    """
    3D epsilon symbol for indices 0, 1, 2.

    epsilon_3d(0,1,2) = +1 (right-handed).
    Returns +1, -1, or 0.
    """
    if (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        return 1
    elif (i, j, k) in [(2, 1, 0), (0, 2, 1), (1, 0, 2)]:
        return -1
    else:
        return 0


def epsilon_nd(indices: Tuple[int, ...]) -> int:
    """
    N-dimensional epsilon symbol. Signature-independent.
    Returns +1, -1, or 0.
    """
    n = len(indices)
    if len(set(indices)) != n:
        return 0
    inversions = 0
    for i in range(n):
        for j in range(i + 1, n):
            if indices[i] > indices[j]:
                inversions += 1
    return 1 if inversions % 2 == 0 else -1


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------

def metric_signature_sign(metric) -> int:
    """
    Return sign(det(metric)) as integer.
    For diag(-,+,+,+): -1
    For diag(+,+,+,+): +1
    """
    try:
        # sympy Matrix
        from sympy import sign as sp_sign
        d = metric.det()
        s = int(sp_sign(d))
        return s
    except AttributeError:
        import numpy as np
        d = float(np.linalg.det(np.array(metric, dtype=float)))
        if d > 0:
            return 1
        elif d < 0:
            return -1
        else:
            return 0


def raise_epsilon_index_sign(metric_inv, indices: tuple) -> int:
    """
    For diagonal orthonormal metric_inv, return product of metric_inv[i,i]
    for i in indices as integer.
    Lorentzian eta_inv, indices (0,1,2,3): (-1)(+1)(+1)(+1) = -1.
    Euclidean delta, indices (0,1,2,3): (+1)(+1)(+1)(+1) = +1.
    """
    result = 1
    for i in indices:
        result *= int(metric_inv[i, i])
    return result


# ---------------------------------------------------------------------------
# Lorentzian epsilon tensor (Phase 1B: orthonormal frame priority)
# ---------------------------------------------------------------------------

def epsilon_tensor_down(metric, a: int, b: int, c: int, d: int) -> int:
    """
    Lower-index epsilon tensor epsilon_abcd.

    For orthonormal frame (|det(g)| = 1):
        epsilon_tensor_down = epsilon_symbol

    For Lorentzian orthonormal frame diag(-,+,+,+):
        epsilon_tensor_down(eta, 0,1,2,3) = +1

    For Euclidean orthonormal frame diag(+,+,+,+):
        epsilon_tensor_down(delta, 0,1,2,3) = +1

    Phase 1B implementation: orthonormal frame (sqrt(|det(g)|) = 1).
    """
    # For orthonormal frames, sqrt(|det(g)|) = 1, so epsilon_tensor_down = epsilon_symbol.
    # metric argument retained for API consistency and future general-metric extension.
    return epsilon_symbol(a, b, c, d)


def epsilon_tensor_up(metric_inv, a: int, b: int, c: int, d: int) -> int:
    """
    Upper-index epsilon tensor epsilon^abcd.

    epsilon^abcd = g^ae g^bf g^cg g^dh epsilon_efgh

    For diagonal orthonormal metric_inv:
        epsilon^abcd = metric_inv[a,a] * metric_inv[b,b] * metric_inv[c,c] * metric_inv[d,d]
                       * epsilon_symbol(a,b,c,d)

    For Lorentzian orthonormal frame diag(-,+,+,+):
        epsilon^0123 = (-1)(+1)(+1)(+1) * (+1) = -1

    For Euclidean orthonormal frame diag(+,+,+,+):
        epsilon^0123 = (+1)(+1)(+1)(+1) * (+1) = +1

    Phase 1B implementation: diagonal orthonormal metric_inv.
    """
    sym = epsilon_symbol(a, b, c, d)
    if sym == 0:
        return 0
    sign = raise_epsilon_index_sign(metric_inv, (a, b, c, d))
    return sign * sym
