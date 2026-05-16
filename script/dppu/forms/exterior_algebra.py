"""
Exterior Differential Form Algebra
===================================

Pure exterior algebra operations on differential forms.

Forms are represented as dictionaries mapping ordered index tuples
to SymPy coefficients::

    Form = Dict[Tuple[int, ...], SymPy_expr]

The empty tuple ``()`` denotes a scalar 0-form.
Index tuples are always kept in ascending sorted order; the
appropriate sign is absorbed into the coefficient when necessary.

Example::

    from dppu.forms.exterior_algebra import basis, wedge, form_to_str
    e0 = basis(0)   # {(0,): 1}
    e1 = basis(1)   # {(1,): 1}
    e01 = wedge(e0, e1)  # {(0, 1): 1}

Author: Muacca
"""

from __future__ import annotations

from typing import Dict, Iterable, Tuple

from sympy import S, factor, simplify

# ---------------------------------------------------------------------------
# Type alias
# ---------------------------------------------------------------------------

Form = Dict[Tuple[int, ...], object]


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def clean_form(form: Form) -> Form:
    """Remove zero-coefficient terms from *form*."""
    out: Form = {}
    for key, val in form.items():
        sval = simplify(val)
        if sval != S.Zero:
            out[key] = sval
    return out


# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

def basis(i: int) -> Form:
    """Return the elementary 1-form ``e^i``."""
    return {(i,): S.One}


def scalar_form(value) -> Form:
    """Wrap *value* as a 0-form.  Returns ``{}`` if *value* is zero."""
    value = simplify(value)
    return {} if value == S.Zero else {(): value}


# ---------------------------------------------------------------------------
# Arithmetic
# ---------------------------------------------------------------------------

def add_forms(*forms: Form) -> Form:
    """Return the sum of any number of forms."""
    out: Form = {}
    for form in forms:
        for key, val in form.items():
            out[key] = out.get(key, S.Zero) + val
    return clean_form(out)


def scale_form(value, form: Form) -> Form:
    """Return *value* · *form*."""
    return clean_form({key: value * val for key, val in form.items()})


def wedge(left: Form, right: Form) -> Form:
    """Return the exterior product ``left ∧ right``."""
    out: Form = {}
    for left_idx, left_val in left.items():
        for right_idx, right_val in right.items():
            if set(left_idx).intersection(right_idx):
                continue  # Repeated index → zero
            seq = list(left_idx) + list(right_idx)
            inversions = sum(
                1
                for i in range(len(seq))
                for j in range(i + 1, len(seq))
                if seq[i] > seq[j]
            )
            sign = -S.One if inversions % 2 else S.One
            key = tuple(sorted(seq))
            out[key] = out.get(key, S.Zero) + sign * left_val * right_val
    return clean_form(out)


# ---------------------------------------------------------------------------
# Display and extraction
# ---------------------------------------------------------------------------

def form_to_str(form: Form) -> str:
    """Human-readable representation of *form*."""
    if not form:
        return "0"
    terms = []
    for key in sorted(form):
        basis_name = "1" if key == () else "e" + "".join(str(i) for i in key)
        terms.append(f"{factor(form[key])}*{basis_name}")
    return " + ".join(terms)


def coefficient(form: Form, key: Tuple[int, ...] = (0, 1, 2, 3)) -> object:
    """Extract and simplify the coefficient for *key* in *form*."""
    return simplify(form.get(key, S.Zero))
