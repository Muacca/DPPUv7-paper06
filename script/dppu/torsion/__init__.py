"""
DPPU Torsion Layer
==================

Provides torsion-related components for Einstein-Cartan gravity.

Modules:
- mode: Torsion ansatz mode enum (AX/VT/MX)
- nieh_yan: Nieh-Yan variant enum and density computation
- ansatz: Torsion tensor construction
- scalar: Torsion scalar T = T_{abc} T^{abc}
"""

from .mode import Mode
from .nieh_yan import NyVariant

__all__ = [
    'Mode',
    'NyVariant',
]
