"""
DPPU Action Layer
=================

Lagrangian, effective potential, EC action builder, and reduced-sector
helpers.  Heavy reduced-sector helpers live in ``dppu.action.reduced_sector``
and are imported explicitly by callers to avoid pipeline import cycles.
"""

from .lagrangian import compute_lagrangian, compute_action
from .potential import compute_effective_potential, get_potential_function, subs_zero_modes
from .ec_action import compute_c2_ec, build_veff_ec, build_veff_ec_func

__all__ = [
    'compute_lagrangian',
    'compute_action',
    'compute_effective_potential',
    'get_potential_function',
    'subs_zero_modes',
    'compute_c2_ec',
    'build_veff_ec',
    'build_veff_ec_func',
]
