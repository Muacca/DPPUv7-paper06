"""
DPPU Forms Package
==================

Exterior differential form algebra and Nieh-Yan form engine for
DPPU-LZ computations.

Modules:
- exterior_algebra: Pure Form type and operations (wedge, add, scale, …)
- nieh_yan:         NiehYanFormEngine — NY form suite for LZ-native S3
"""

from .exterior_algebra import (
    Form,
    basis,
    scalar_form,
    clean_form,
    add_forms,
    scale_form,
    wedge,
    form_to_str,
    coefficient,
)
from .nieh_yan import NiehYanFormEngine

__all__ = [
    # exterior_algebra
    "Form",
    "basis",
    "scalar_form",
    "clean_form",
    "add_forms",
    "scale_form",
    "wedge",
    "form_to_str",
    "coefficient",
    # nieh_yan
    "NiehYanFormEngine",
]
