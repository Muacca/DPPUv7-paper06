"""
Canonical Spin-2 Basis Utilities
================================

Provides the canonical traceless symmetric spin-2 basis matrices T1-T5
used for metric perturbations, and associated validation routines.

Author: Muacca
"""

from typing import Dict, Optional

from sympy import Matrix, S, simplify, sqrt


def build_spin2_basis() -> Dict[str, Matrix]:
    """
    Return the canonical traceless symmetric spin-2 basis T1..T5.
    
    The basis is orthonormal under the trace inner product: tr(T_A T_B) = \delta_{AB}
    and traceless: tr(T_A) = 0.
    """
    return {
        "T1": Matrix.diag(1, 1, -2) / sqrt(6),
        "T2": Matrix.diag(1, -1, 0) / sqrt(2),
        "T3": Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / sqrt(2),
        "T4": Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / sqrt(2),
        "T5": Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / sqrt(2),
    }


def validate_spin2_basis(basis: Optional[Dict[str, Matrix]] = None) -> bool:
    """
    Assert tr(T_A)=0 and tr(T_A T_B)=delta_AB.
    
    Raises AssertionError if the basis is invalid.
    """
    basis = build_spin2_basis() if basis is None else basis
    labels = ["T1", "T2", "T3", "T4", "T5"]
    
    for label in labels:
        tr_val = simplify(basis[label].trace())
        if tr_val != S.Zero:
            raise AssertionError(f"{label} trace is not zero: {tr_val}")
            
    for a_label in labels:
        for b_label in labels:
            inner = simplify((basis[a_label] * basis[b_label]).trace())
            expected = S.One if a_label == b_label else S.Zero
            if inner != expected:
                raise AssertionError(
                    f"{a_label},{b_label} inner product mismatch: {inner}"
                )
    return True


def beta_squash_matrix() -> Matrix:
    """
    Return the unnormalized beta squash generator.
    
    This corresponds to the diag(1, 1, -2) deformation.
    """
    return Matrix.diag(1, 1, -2)
