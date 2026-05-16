"""
Symbolic Computation Utilities
==============================

Provides utilities for symbolic zero-proving and numerical witness finding.

Used primarily for Riemann tensor antisymmetry verification (3-level approach):
1. Try to prove symbolically that an expression is zero
2. If not proved, try to find a numerical witness that it's nonzero
3. Based on results, classify as PROVED_ZERO, WITNESS_NONZERO, or UNPROVED

Author: Muacca
"""

from typing import Dict, List, Optional, Tuple, Any

import numpy as np
from sympy import (
    simplify, expand, factor, cancel, together, trigsimp, Symbol
)
from sympy.core.expr import Expr


def prove_zero(
    expr: Expr,
    assumptions_dict: Optional[Dict] = None,
    timeout_seconds: int = 10
) -> Tuple[bool, str]:
    """
    Try to prove that a symbolic expression equals zero.

    Uses a multi-stage simplification pipeline:
    1. Lightweight: simplify, expand
    2. Heavy: factor, cancel, together, trigsimp
    3. Combined: pipeline of all

    Args:
        expr: SymPy expression to test
        assumptions_dict: Optional substitutions to apply first
        timeout_seconds: (Currently unused) Timeout for heavy operations

    Returns:
        Tuple of (proved: bool, method: str)
        - If proved: (True, "PROVED_ZERO (via <method>)")
        - If not: (False, "UNPROVED")

    Examples:
        >>> from sympy import symbols, sin, cos
        >>> x = symbols('x')
        >>> prove_zero(sin(x)**2 + cos(x)**2 - 1)
        (True, 'PROVED_ZERO (via trigsimp)')
    """
    # Trivial case
    if expr == 0:
        return True, "PROVED_ZERO (trivial)"

    # Apply substitutions if provided
    if assumptions_dict:
        expr = expr.subs(assumptions_dict)

    # Stage 0: cancel() — immediate decision for rational functions (z_i = exp(q_i/√2) rationalized)
    # For expressions without transcendental functions, cancel() acts as a complete algorithm,
    # allowing the subsequent expensive simplify/trigsimp steps to be skipped.
    try:
        if cancel(expr) == 0:
            return True, "PROVED_ZERO (via cancel)"
    except Exception:
        pass

    # Stage 1: Lightweight simplifications
    for name, func in [("simplify", simplify), ("expand", expand)]:
        try:
            result = func(expr)
            if result == 0:
                return True, f"PROVED_ZERO (via {name})"
        except Exception:
            continue

    # Stage 2: Heavy simplifications
    for name, func in [
        ("factor", factor),
        ("cancel", cancel),
        ("together", together),
        ("trigsimp", trigsimp)
    ]:
        try:
            result = func(expr)
            if result == 0:
                return True, f"PROVED_ZERO (via {name})"
        except Exception:
            continue

    # Stage 3: Combined pipeline
    try:
        result = trigsimp(simplify(factor(together(expr))))
        if result == 0:
            return True, "PROVED_ZERO (via pipeline)"
    except Exception:
        pass

    return False, "UNPROVED"


def generate_test_points(
    symbols_list: List[Symbol],
    n_points: int = 10
) -> List[Dict[Symbol, float]]:
    """
    Generate diverse test points for numerical evaluation.

    Avoids special values (0, 1, -1) to reduce false negatives.

    Args:
        symbols_list: List of SymPy symbols to generate values for
        n_points: Number of test points to generate

    Returns:
        List of dicts mapping symbols to float values

    Notes:
        - Positive symbols get values in (0.1, 10.0)
        - Other symbols get values in (-10.0, 10.0)
    """
    points = []
    for _ in range(n_points):
        point = {}
        for sym in symbols_list:
            if hasattr(sym, 'is_positive') and sym.is_positive:
                point[sym] = np.random.uniform(0.1, 10.0)
            else:
                point[sym] = np.random.uniform(-10.0, 10.0)
        points.append(point)
    return points


def find_nonzero_witness(
    expr: Expr,
    symbols_list: List[Symbol],
    n_points: int = 10,
    precision: int = 50
) -> Tuple[bool, Optional[Dict], Optional[float]]:
    """
    Try to find a numerical witness that an expression is nonzero.

    Uses high-precision arithmetic (mpmath) to evaluate the expression
    at random test points.

    Args:
        expr: SymPy expression to test
        symbols_list: List of symbols in the expression
        n_points: Number of test points to try
        precision: Number of decimal places for mpmath

    Returns:
        Tuple of (found: bool, witness_point: dict or None, witness_value: float or None)
        - If witness found: (True, {sym: val, ...}, value)
        - If all points ~0: (False, None, None)

    Notes:
        - A threshold of 10^(-precision+10) is used to detect nonzero values
        - Singular points (division by zero, etc.) are skipped
    """
    try:
        import mpmath
        from sympy import lambdify
        f = lambdify(symbols_list, expr, modules=['mpmath'])
    except Exception:
        return False, None, None

    test_points = generate_test_points(symbols_list, n_points)
    mpmath.mp.dps = precision
    threshold = mpmath.mpf(10) ** (-precision + 10)

    for point in test_points:
        try:
            val = f(*point.values())
            if abs(val) > threshold:
                return True, point, float(val)  # Witness found
        except (ValueError, ZeroDivisionError, TypeError):
            continue  # Skip singular points

    return False, None, None  # All points ~0


def normalize_expression(expr: Expr) -> Expr:
    """
    Normalize a symbolic expression for comparison.

    Applies a standard sequence of simplifications to bring
    equivalent expressions to a canonical form.

    Args:
        expr: SymPy expression

    Returns:
        Normalized expression

    Notes:
        This is useful for comparing residuals across different
        tensor components.
    """
    try:
        result = simplify(expand(expr))
        result = factor(together(result))
        return result
    except Exception:
        return expr


def derivative_inventory(
    expr: Expr,
    symbol_pairs: "List[Tuple[str, Any]]",
) -> "Dict[str, bool]":
    """
    Check which symbols appear as free symbols in *expr*.

    Args:
        expr:         SymPy expression to inspect.
        symbol_pairs: Iterable of (name, symbol) pairs.  Each *symbol*
                      is tested for membership in ``expr.free_symbols``.

    Returns:
        Dict mapping each *name* to True/False.

    Example::

        from sympy import symbols
        x, y = symbols("x y")
        derivative_inventory(x**2, [("x", x), ("y", y)])
        # → {"x": True, "y": False}
    """
    return {name: (sym in expr.free_symbols) for name, sym in symbol_pairs}


def nonzero_component_count(
    array_expr,
    subs: "Optional[Dict]" = None,
) -> "Tuple[int, List]":
    """
    Count the nonzero components of a symbolic tensor array.

    Rank and dimension are inferred from ``array_expr.shape``.

    Args:
        array_expr: ``MutableDenseNDimArray`` (or any object supporting
                    integer-tuple indexing and ``.shape``).
        subs:       Optional substitution dict applied before testing.

    Returns:
        ``(count, examples)`` where *count* is the number of nonzero
        components and *examples* is a list of up to eight
        ``((index_tuple), factor(value))`` pairs.
    """
    from sympy import simplify as _simplify, factor as _factor

    shape = array_expr.shape
    rank = len(shape)
    ranges = [range(d) for d in shape]

    count = 0
    examples: List = []

    def _visit(prefix: List[int], depth: int) -> None:
        nonlocal count
        if depth == rank:
            val = array_expr[tuple(prefix)]
            if subs:
                val = val.subs(subs)
            val = _simplify(val)
            if val != 0:
                count += 1
                if len(examples) < 8:
                    examples.append((tuple(prefix), _factor(val)))
            return
        for i in ranges[depth]:
            _visit(prefix + [i], depth + 1)

    _visit([], 0)
    return count, examples
