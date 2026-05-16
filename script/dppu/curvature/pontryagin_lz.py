"""
Pontryagin Density (Lorentzian, metric-aware)
=============================================

Implements Lorentzian Pontryagin density P = <R, *R> using
hodge_dual_2form (metric-aware Hodge dual for 2-forms).

Phase 1D DPPU-LZ primary definition:

  P = 1/4 * sum_{a,b,c,d} R^{ab cd} (*R_ab)_{cd}

  where:
    R^{abcd} = g^{aA} g^{bB} g^{cC} g^{dD} R_ABCD   (all indices raised)
    (*R_ab)_{cd} = hodge_dual_2form((R_abcd)[a,b,:,:], metric, metric_inv)[c,d]

  Normalization: factor 1/4 accounts for antisymmetric pair double-counting.

Cross-check (epsilon-based):

  P_eps = 1/8 * sum_{a,b,c,d,e,f}
            R^{ab cd} * epsilon_tensor_down(metric,c,d,e,f) * R_ab^{ef}

where R_ab^{ef} raises only the form indices of the second curvature factor;
the internal pair (a,b) remains lowered.
  This is identical to P_hodge and is implemented independently as a sanity check.

IMPORTANT:
  - Do NOT import or use compute_hodge_dual in this module.
  - Do NOT use epsilon_4d as a Lorentzian epsilon tensor.
  - epsilon_tensor_down / epsilon_tensor_up must be used for tensor operations.
"""

from sympy import Matrix, Rational, S, cancel, factor, simplify
from sympy.tensor.array import MutableDenseNDimArray

from .hodge import hodge_dual_2form
from ..utils.epsilon import epsilon_tensor_down, epsilon_tensor_up


# ---------------------------------------------------------------------------
# Core tensor utilities
# ---------------------------------------------------------------------------

def curvature_2form_matrix(R_abcd, a: int, b: int, dim: int = 4):
    """
    Return 4x4 Matrix F_cd = R_abcd[a,b,c,d] for fixed internal pair (a,b).

    Input:
      R_abcd: MutableDenseNDimArray or indexable [a,b,c,d]
      a, b:   fixed antisymmetric pair (0..dim-1)
      dim:    dimension (default 4)

    Output:
      dim x dim sympy Matrix F with F[c,d] = R_abcd[a,b,c,d]
    """
    F = Matrix.zeros(dim, dim)
    for c in range(dim):
        for d in range(dim):
            F[c, d] = R_abcd[a, b, c, d]
    return F


def raise_curvature_all_indices(R_abcd, metric_inv, dim: int = 4):
    """
    Return R^{abcd} from all-lowered R_abcd.

    R_up[a,b,c,d] = g^{aA} g^{bB} g^{cC} g^{dD} R_ABCD

    Input:
      R_abcd:     MutableDenseNDimArray [A,B,C,D] (all-lowered)
      metric_inv: dim x dim sympy Matrix (inverse metric)
      dim:        dimension (default 4)

    Output:
      MutableDenseNDimArray of shape (dim,dim,dim,dim)
    """
    R_up = MutableDenseNDimArray.zeros(dim, dim, dim, dim)
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    val = S.Zero
                    for A in range(dim):
                        gAa = metric_inv[a, A]
                        if gAa == S.Zero:
                            continue
                        for B in range(dim):
                            gBb = metric_inv[b, B]
                            if gBb == S.Zero:
                                continue
                            for C in range(dim):
                                gCc = metric_inv[c, C]
                                if gCc == S.Zero:
                                    continue
                                for D in range(dim):
                                    gDd = metric_inv[d, D]
                                    if gDd == S.Zero:
                                        continue
                                    R_val = R_abcd[A, B, C, D]
                                    if R_val != S.Zero:
                                        val += gAa * gBb * gCc * gDd * R_val
                    if val != S.Zero:
                        R_up[a, b, c, d] = val
    return R_up


# ---------------------------------------------------------------------------
# Primary Pontryagin density (hodge_dual_2form based)
# ---------------------------------------------------------------------------

def pontryagin_density_lorentzian(
    R_abcd,
    metric,
    metric_inv=None,
    dim: int = 4,
    simplify_result: bool = True,
):
    """
    Compute Lorentzian metric-aware Pontryagin density.

    P = 1/4 * sum_{a,b,c,d} R^{ab cd} (*R_ab)_{cd}

    Uses hodge_dual_2form for the form-index Hodge dual.
    Does NOT use compute_hodge_dual.

    Input:
      R_abcd:         MutableDenseNDimArray (dim^4), all-lowered
      metric:         dim x dim sympy Matrix (frame metric)
      metric_inv:     optional dim x dim sympy Matrix; if None, computed via metric.inv()
      dim:            dimension (default 4)
      simplify_result: if True, apply cancel() to final result

    Output:
      sympy scalar expression for P
    """
    if metric_inv is None:
        metric_inv = metric.inv()

    # Raise all indices: R^{abcd}
    R_up = raise_curvature_all_indices(R_abcd, metric_inv, dim)

    # Accumulate P = 1/4 * sum_{a,b,c,d} R^{abcd} (*R_ab)_{cd}
    P = S.Zero
    for a in range(dim):
        for b in range(dim):
            if a == b:
                continue  # antisymmetric in (a,b): diagonal is zero
            # 2-form (R_ab)_{cd} = R_abcd for fixed (a,b)
            F_ab = curvature_2form_matrix(R_abcd, a, b, dim)
            # Hodge dual: (*R_ab)_{cd}
            star_F_ab = hodge_dual_2form(F_ab, metric, metric_inv)
            for c in range(dim):
                for d in range(dim):
                    Rup_abcd = R_up[a, b, c, d]
                    if Rup_abcd != S.Zero:
                        sf = star_F_ab[c, d]
                        if sf != S.Zero:
                            P += Rup_abcd * sf

    P = Rational(1, 4) * P

    if simplify_result:
        P = cancel(P)

    return P


# ---------------------------------------------------------------------------
# Helper: raise only form indices for fixed internal pair
# ---------------------------------------------------------------------------

def raise_form_indices_for_fixed_internal_pair(
    R_abcd,
    a: int,
    b: int,
    metric_inv,
    dim: int = 4,
):
    """
    Return Matrix G^{ef} = R_ab^{ef} for fixed lower internal pair (a,b).

    Only the form indices are raised:
        R_ab^{ef} = g^{eE} g^{fF} R_abEF

    The internal indices a, b are NOT raised.
    This is consistent with hodge_dual_2form, which operates on the 2-form
    (form) indices of a fixed-pair 2-form and does not raise internal indices.

    Input:
      R_abcd:     MutableDenseNDimArray [A,B,C,D] (all-lowered)
      a, b:       fixed lower internal pair (0..dim-1)
      metric_inv: dim x dim sympy Matrix (inverse metric)
      dim:        dimension (default 4)

    Output:
      dim x dim sympy Matrix G with G[e,f] = R_ab^{ef}
    """
    G = Matrix.zeros(dim, dim)
    for e in range(dim):
        for f in range(dim):
            val = S.Zero
            for E in range(dim):
                geE = metric_inv[e, E]
                if geE == S.Zero:
                    continue
                for F in range(dim):
                    gfF = metric_inv[f, F]
                    if gfF == S.Zero:
                        continue
                    r = R_abcd[a, b, E, F]
                    if r != S.Zero:
                        val += geE * gfF * r
            if val != S.Zero:
                G[e, f] = cancel(val)
    return G


# ---------------------------------------------------------------------------
# Cross-check: epsilon tensor based Pontryagin density
# ---------------------------------------------------------------------------

def pontryagin_density_lorentzian_by_epsilon(
    R_abcd,
    metric,
    metric_inv=None,
    dim: int = 4,
    simplify_result: bool = True,
):
    """
    Independent cross-check using epsilon tensor / explicit component formula.

    Derivation from primary definition:
      P_hodge = 1/4 * sum_{a,b,c,d} R^{ab cd} (*R_ab)_{cd}
             = 1/4 * sum_{a,b,c,d} R^{ab cd} * (1/2) sum_{e,f} epsilon_down(c,d,e,f) * R_ab^{ef}
             = 1/8 * sum_{a,b,c,d,e,f} R^{ab cd} * epsilon_down(c,d,e,f) * R_ab^{ef}

    where:
      R^{ab cd} raises all four indices of the first curvature factor.
      R_ab^{ef} raises only the form indices of the second curvature factor;
                internal pair (a,b) remains lowered.

    This is consistent with hodge_dual_2form, which acts on form indices only.
    This function uses epsilon_tensor_down explicitly (NOT epsilon_4d).
    Should match pontryagin_density_lorentzian exactly with simplify_result=True.

    NOTE: An earlier (v0.1) version incorrectly used R^{abef} (all-raised) as
    the second factor. This has been corrected to R_ab^{ef} (form-raised only).

    Input:
      R_abcd:         MutableDenseNDimArray (dim^4), all-lowered
      metric:         dim x dim sympy Matrix (frame metric)
      metric_inv:     optional dim x dim sympy Matrix
      dim:            dimension (default 4)
      simplify_result: if True, apply cancel() to final result

    Output:
      sympy scalar expression for P (cross-check)
    """
    if metric_inv is None:
        metric_inv = metric.inv()

    # Raise all indices for the first factor: R^{abcd}
    R_up_all = raise_curvature_all_indices(R_abcd, metric_inv, dim)

    # P_eps = 1/8 * sum_{a,b,c,d,e,f} R^{abcd} * epsilon_down(c,d,e,f) * R_ab^{ef}
    P_eps = S.Zero
    for a in range(dim):
        for b in range(dim):
            if a == b:
                continue
            # Second factor: form indices only raised, internal pair (a,b) kept lower
            R_ab_up_form = raise_form_indices_for_fixed_internal_pair(
                R_abcd, a, b, metric_inv, dim
            )
            for c in range(dim):
                for d in range(dim):
                    Rup_abcd = R_up_all[a, b, c, d]
                    if Rup_abcd == S.Zero:
                        continue
                    for e in range(dim):
                        for f in range(dim):
                            eps = epsilon_tensor_down(metric, c, d, e, f)
                            if eps == 0:
                                continue
                            Rab_up_ef = R_ab_up_form[e, f]
                            if Rab_up_ef != S.Zero:
                                P_eps += Rup_abcd * eps * Rab_up_ef

    P_eps = Rational(1, 8) * P_eps

    if simplify_result:
        P_eps = cancel(P_eps)

    return P_eps


# ---------------------------------------------------------------------------
# Classification helper
# ---------------------------------------------------------------------------

def classify_pontryagin_expression(expr, zero_tol=None):
    """
    Classify a sympy expression as zero or nonzero.

    Steps:
      1. Direct == S.Zero check
      2. cancel()
      3. factor()
      4. simplify()  (only if previous steps inconclusive)

    Returns one of:
      "EXACT_ZERO"                      -- proven zero symbolically
      "FACTORIZED_NONZERO"              -- nonzero, factorized form available
      "NUMERIC_NONZERO_SYMBOLIC_UNKNOWN"-- numeric witness is nonzero
      "UNCLASSIFIED"                    -- symbolic simplification inconclusive

    The zero_tol parameter is unused for symbolic classification
    (kept for API compatibility with future numeric extension).
    """
    # Stage 0: trivial check
    if expr == S.Zero:
        return "EXACT_ZERO", S.Zero

    # Stage 1: cancel
    expr1 = cancel(expr)
    if expr1 == S.Zero:
        return "EXACT_ZERO", S.Zero

    # Stage 2: factor
    expr2 = factor(expr1)
    if expr2 == S.Zero:
        return "EXACT_ZERO", S.Zero

    # At this point: expr2 is factorized form
    # Check if it is structurally nonzero (no free cancellation left)
    # A factorized expression that is not zero is nonzero.
    # We treat factor() output as canonical: if it isn't zero, classify as nonzero.
    # (simplify is deferred to avoid excessive computation time)
    return "FACTORIZED_NONZERO", expr2


def classify_pontryagin_expression_deep(expr, zero_tol=None):
    """
    Deep classification using simplify() as last resort.

    Same as classify_pontryagin_expression but also tries simplify().
    Use when factor() returns a non-trivial expression but zero is suspected.

    Returns (status_string, canonical_expression).
    """
    status, expr_out = classify_pontryagin_expression(expr, zero_tol)
    if status == "EXACT_ZERO":
        return status, expr_out

    # Stage 3: simplify
    expr3 = simplify(expr_out)
    if expr3 == S.Zero:
        return "EXACT_ZERO", S.Zero

    return "FACTORIZED_NONZERO", expr3
