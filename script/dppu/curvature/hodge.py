"""
Hodge Dual Operator
===================

Provides metric-aware Hodge dual for 2-forms (DPPU-LZ native convention).

Convention (DPPU-LZ v1.1, Lorentzian signature (-,+,+,+)):
  (*F)_ab = (1/2) epsilon_tensor_down(metric, a,b,c,d) * F^cd
  where F^cd = metric_inv^{ce} metric_inv^{df} F_ef

  Hodge basis table (Lorentzian, epsilon_0123 = +1):
    *(01) = -(23)
    *(02) = +(13)
    *(03) = -(12)
    *(12) = +(03)
    *(13) = -(02)
    *(23) = +(01)

  Double Hodge on 2-forms: ** = -1  (Lorentzian)

IMPORTANT: This module does NOT import legacy_index.
           Legacy index comparisons belong in check scripts only.
"""

import numpy as np
from ..utils.epsilon import epsilon_symbol, epsilon_tensor_down


def hodge_dual_2form(F_down, metric, metric_inv=None, signature="lorentzian"):
    """
    Compute lower-index Hodge dual of a 2-form.

    (*F)_ab = (1/2) sum_{c,d} epsilon_tensor_down(metric, a,b,c,d) * F^cd

    where F^cd = sum_{e,f} metric_inv^{ce} metric_inv^{df} F_ef

    Input:
      F_down: antisymmetric 4x4 sympy Matrix with F_down[a,b] = F_ab
      metric: 4x4 sympy Matrix (frame metric)
      metric_inv: optional 4x4 sympy Matrix (inverse metric);
                  if None, computed via metric.inv()
      signature: "lorentzian" or "euclidean" (for logging/assertion only;
                 actual calculation uses metric)

    Output:
      star_F_down: 4x4 sympy Matrix with (star_F)_ab
    """
    from sympy import Matrix, Rational

    if metric_inv is None:
        metric_inv = metric.inv()

    # Raise indices: F^cd = sum_{e,f} metric_inv[c,e] * metric_inv[d,f] * F_down[e,f]
    F_up = [[0] * 4 for _ in range(4)]
    for c in range(4):
        for d in range(4):
            val = 0
            for e in range(4):
                for f in range(4):
                    mce = metric_inv[c, e]
                    mdf = metric_inv[d, f]
                    if mce != 0 and mdf != 0:
                        val += mce * mdf * F_down[e, f]
            F_up[c][d] = val

    # (*F)_ab = (1/2) sum_{c,d} epsilon_tensor_down(metric, a,b,c,d) * F^cd
    star_F = [[0] * 4 for _ in range(4)]
    for a in range(4):
        for b in range(4):
            val = 0
            for c in range(4):
                for d in range(4):
                    eps = epsilon_tensor_down(metric, a, b, c, d)
                    if eps != 0 and F_up[c][d] != 0:
                        val += Rational(1, 2) * eps * F_up[c][d]
            star_F[a][b] = val

    return Matrix(star_F)


def basis_2form(a: int, b: int):
    """
    Return antisymmetric 4x4 sympy Matrix for basis 2-form e^a wedge e^b.
    Components: F[a,b] = +1, F[b,a] = -1, all others 0.
    """
    from sympy import Matrix
    F = Matrix.zeros(4, 4)
    F[a, b] = 1
    F[b, a] = -1
    return F


# ---------------------------------------------------------------------------
# Legacy API (kept for backward compatibility; uses epsilon_symbol as symbol)
# ---------------------------------------------------------------------------

def compute_hodge_dual(R_tensor: np.ndarray) -> np.ndarray:
    """
    [DEPRECATED] Compute Hodge dual on (cd) indices of R^{ab}_{cd}.

    (*R)^{ab}_{cd} = (1/2) epsilon_symbol(c,d,e,f) R^{ab,ef}

    WARNING: Uses epsilon_symbol (combinatorial, NOT metric-aware tensor).
    This function does NOT account for Lorentzian metric signature.
    For Lorentzian / metric-aware Hodge on 2-forms, use hodge_dual_2form().
    Using this function for Pontryagin density evaluation will give WRONG signs.

    Args:
        R_tensor: numpy array (4,4,4,4) representing R^{ab}_{cd}

    Returns:
        numpy array (4,4,4,4) representing (*R)^{ab}_{cd}
    """
    import warnings
    warnings.warn(
        "compute_hodge_dual is non-metric-aware (uses epsilon_symbol only). "
        "Use hodge_dual_2form() for Lorentzian computations.",
        DeprecationWarning,
        stacklevel=2,
    )
    R_dual = np.zeros_like(R_tensor)

    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    val = 0.0
                    for e in range(4):
                        for f in range(4):
                            eps = epsilon_symbol(c, d, e, f)
                            if eps != 0:
                                val += 0.5 * eps * R_tensor[a, b, e, f]
                    R_dual[a, b, c, d] = val

    return R_dual


def classify_lz_2form_block(c: int, d: int) -> str:
    """
    DPPU-LZ native 2-form block classification.
    time-space: one index is 0 and the other is in {1,2,3}
    spatial: both indices are in {1,2,3}
    zero: c == d (diagonal, antisymmetric form is zero)
    """
    if c == d:
        return "zero"
    if 0 in (c, d):
        return "time-space"
    return "spatial"


def classify_block_legacy_euclidean(c: int, d: int) -> str:
    """
    [DEPRECATED] Legacy Euclidean DPPU block classification.
    Spatial: c,d in {0,1,2}; fiber/mixed: involves index 3.
    Do NOT use for DPPU-LZ native Lorentzian computations.
    Use classify_lz_2form_block() instead.
    """
    if c < 3 and d < 3:
        return 'spatial'
    return 'mixed'


def hodge_swaps_blocks() -> bool:
    """
    Property: Lorentzian Hodge dual swaps time-space and spatial blocks.
    * : time-space -> spatial
    * : spatial -> time-space
    """
    return True
