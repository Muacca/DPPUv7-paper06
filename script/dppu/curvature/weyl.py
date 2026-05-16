"""
Weyl Tensor and Scalar Computation
==================================

Computes the Weyl tensor (conformal curvature) and Weyl scalar (C^2).

Mathematical Definition:
    C_{abcd} = R_{abcd} 
             - (1/2)(g_{ac}R_{bd} - g_{ad}R_{bc} - g_{bc}R_{ad} + g_{bd}R_{ac})
             + (R/6)(g_{ac}g_{bd} - g_{ad}g_{bc})

    C^2 = C_{abcd}C^{abcd}

Properties:
    - Trace-free: C^a_{bad} = 0
    - Same symmetries as Riemann tensor
    - Vanishes for conformally flat spacetimes (e.g., S^3 x R is conformally flat?)
      Actually S^3 is conformally flat, but S^3 x S^1 is not necessarily so.
      However, S^3 x R^1 is conformally flat.

Author: Muacca
"""

from typing import Any, Optional, Tuple

from sympy import S, cancel, Matrix, Rational
from sympy.tensor.array import MutableDenseNDimArray


def compute_weyl_tensor(
    R_abcd: MutableDenseNDimArray,
    Ricci: Matrix,
    R_scalar: Any,
    metric: Matrix,
    dim: int,
    logger: Optional[Any] = None
) -> MutableDenseNDimArray:
    """
    Compute Weyl tensor C_{abcd}.

    Args:
        R_abcd: Riemann tensor with all lower indices
        Ricci: Ricci tensor R_{ab}
        R_scalar: Ricci scalar R
        metric: Metric tensor g_{ab} (or eta_{ab})
        dim: Dimension (must be >= 3, usually 4)
        logger: Optional logger

    Returns:
        Weyl tensor C_{abcd} as MutableDenseNDimArray
    """
    if logger:
        logger.info("Computing Weyl Tensor C_abcd...")

    if dim < 3:
        # Weyl tensor vanishes for dim < 3
        return MutableDenseNDimArray.zeros(dim, dim, dim, dim)

    C_abcd = MutableDenseNDimArray.zeros(dim, dim, dim, dim)

    # Precompute R/((d-1)(d-2)) factor?
    # Standard 4D formula:
    # C = R - 1/2 (gR + gR) + R/6 (gg)

    # Generalized formula for dimension n:
    # C_{abcd} = R_{abcd} 
    #          - (1/(n-2)) * (g_{ac}R_{bd} - g_{ad}R_{bc} - g_{bc}R_{ad} + g_{bd}R_{ac})
    #          + (R/((n-1)(n-2))) * (g_{ac}g_{bd} - g_{ad}g_{bc})
    
    # Using n=4 for this implementation as per requirement
    n = dim
    if n == 4:
        coeff_Ricci = Rational(1, 2)
        coeff_Scalar = Rational(1, 6)
    else:
        coeff_Ricci = Rational(1, n - 2)
        coeff_Scalar = Rational(1, (n - 1) * (n - 2))

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    # Riemann term
                    term = R_abcd[a, b, c, d]
                    
                    # Ricci term: - 1/2 (g_{ac}R_{bd} - g_{ad}R_{bc} - g_{bc}R_{ad} + g_{bd}R_{ac})
                    # Note signs:
                    # Let's follow the standard definition carefully.
                    # C_{abcd} = R_{abcd} + (1/(n-2)) * (g_{ad}R_{bc} + g_{bc}R_{ad} - g_{ac}R_{bd} - g_{bd}R_{ac})
                    # ... wait, let's verify standard formula signs.
                    # R_{abcd} is antisymmetric in [ab] and [cd].
                    # The term g_{ac}R_{bd} is not antisymmetric in [ab]. 
                    # Correct structure is (g_{ac}R_{bd} - g_{ad}R_{bc} + g_{bd}R_{ac} - g_{bc}R_{ad}) ???
                    # Let's check Wikipedia or standard text.
                    # Wikipedia:
                    # C_{iklm} = R_{iklm} + (1/(n-2))(R_{im}g_{kl} - R_{il}g_{km} + R_{kl}g_{im} - R_{km}g_{il}) ...
                    # This depends on sign convention of Riemann.
                    
                    # My Riemann definition: R^a_{bcd} = ...
                    # R_{abcd} = g_{ae} R^e_{bcd}
                    
                    # Let's use the explicit formula from plan:
                    # C_{abcd} = R_{abcd} - 1/2(g_{ac}R_{bd} - g_{ad}R_{bc} - g_{bc}R_{ad} + g_{bd}R_{ac})
                    #            + R/6(g_{ac}g_{bd} - g_{ad}g_{bc})
                    
                    term_Ricci = (
                        metric[a, c] * Ricci[b, d] -
                        metric[a, d] * Ricci[b, c] -
                        metric[b, c] * Ricci[a, d] +
                        metric[b, d] * Ricci[a, c]
                    )
                    
                    term -= coeff_Ricci * term_Ricci
                    
                    # Scalar term
                    term_Scalar = (
                        metric[a, c] * metric[b, d] -
                        metric[a, d] * metric[b, c]
                    )
                    
                    term += coeff_Scalar * R_scalar * term_Scalar
                    
                    if term != 0:
                        C_abcd[a, b, c, d] = cancel(term)

    if logger:
        logger.success("Weyl Tensor computed")

    return C_abcd


def compute_weyl_scalar(
    C_abcd: MutableDenseNDimArray,
    metric_inv: Matrix,
    dim: int,
    logger: Optional[Any] = None
) -> Any:
    """
    Compute Weyl scalar C^2 = C_{abcd}C^{abcd}.

    Args:
        C_abcd: Weyl tensor with all lower indices
        metric_inv: Inverse metric g^{ab} (or eta^{ab})
        dim: Dimension
        logger: Optional logger

    Returns:
        Weyl scalar C^2 as SymPy expression
    """
    if logger:
        logger.info("Computing Weyl Scalar C^2...")

    # C^2 = g^{ae} g^{bf} g^{cg} g^{dh} C_{abcd} C_{efgh}
    
    # Optimization: Since C_{abcd} is sparse/symbolic, we can iterate non-zero elements?
    # Or just iterate all indices. 4^4 = 256 iterations, acceptable.
    
    # First, raise indices for one copy of C
    # C^{abcd} = g^{ae} g^{bf} g^{cg} g^{dh} C_{efgh}
    
    # Actually, let's do C_{abcd} C^{abcd} directly by summation
    
    C_sq = S.Zero
    
    # To reduce nesting depth, maybe pre-calculate raised tensor?
    # But C_abcd is easier to pass.
    
    # 4 loops for first tensor, 4 loops for second tensor is too much (4^8).
    # We must contract carefully.
    
    # C^2 = C_{abcd} g^{ae} g^{bf} g^{cg} g^{dh} C_{efgh}
    
    # Let's define C_raised
    C_upper = MutableDenseNDimArray.zeros(dim, dim, dim, dim)
    
    has_nonzero = False
    
    # Raise one index at a time? No, just do it.
    # C^{abcd}
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    val = S.Zero
                    # contraction C^{abcd} = sum_{efgh} g^{ae} g^{bf} g^{cg} g^{dh} C_{efgh}
                    # This is still 4^4 loops inside 4^4 loops.
                    
                    # Optimization for diagonal metric (frame metric)
                    # eta = diag(1, 1, 1, 1) usually.
                    # if metric is diagonal:
                    # C^{abcd} = g^{aa} g^{bb} g^{cc} g^{dd} C_{abcd} (no sum)
                    
                    # Let's assume diagonal metric for performance, but check if it is diagonal.
                    is_diagonal = True
                    for i in range(dim):
                        for j in range(dim):
                            if i != j and metric_inv[i, j] != 0:
                                is_diagonal = False
                                break
                    
                    if is_diagonal:
                         # For diagonal metric (frame metric):
                         # C^{abcd} = inv[a,a]*inv[b,b]*inv[c,c]*inv[d,d] * C_{abcd}
                         val = (metric_inv[a, a] * metric_inv[b, b] * 
                                metric_inv[c, c] * metric_inv[d, d] * 
                                C_abcd[a, b, c, d])
                    else:
                        # Full contraction fallback
                        for e in range(dim):
                            for f in range(dim):
                                for g in range(dim):
                                    for h in range(dim):
                                        term = (metric_inv[a, e] * metric_inv[b, f] * 
                                                metric_inv[c, g] * metric_inv[d, h] * 
                                                C_abcd[e, f, g, h])
                                        val += term
                    
                    if val != 0:
                        C_upper[a, b, c, d] = cancel(val)
                        has_nonzero = True

    # Now contract C_{abcd} C^{abcd}
    if not has_nonzero:
        return S.Zero

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    term = C_abcd[a, b, c, d] * C_upper[a, b, c, d]
                    C_sq += term

    C_sq = cancel(C_sq)

    if logger:
        logger.info(f"  Weyl Scalar C^2 = {C_sq}")

    return C_sq
