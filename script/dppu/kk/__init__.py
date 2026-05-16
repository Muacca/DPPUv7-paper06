"""
DPPU KK: Kaluza-Klein Photon Effective Theory
==============================================

Provides the two-route pipeline for computing the O(A²) photon effective
action from 5D Einstein-Cartan gravity.

Routes
------
Shortcut route (fast):
    field_strength → gamma_gamma → extractor

    1. Build F̃_{jk} (topology-aware modified field strength)
    2. Build ω⁽¹⁾ (O(A¹) connection perturbation via Koszul)
    3. Compute R^{Γ²} = Σ ω⁽¹⁾×ω⁽¹⁾  (Γ×Γ shortcut)
    4. Extract Maxwell / CS / mass² coefficients

Validation route (slow):
    field_strength + full_riemann + extractor

    1-2. Same as above
    3b. Compute full R[ω⁽⁰⁾+ω⁽¹⁾] and extract O(A²) piece
    4.  Compare with shortcut result → should be identical

Usage
-----
Quick start::

    from sympy import symbols, Symbol
    from dppu.kk import (
        s3_corrections, make_F_tilde, make_omega1,
        gamma_gamma_ricci, extract_all,
    )

    r0, L = symbols('r0 L', positive=True)
    A   = [Symbol('A0'), Symbol('A1'), Symbol('A2')]
    dA  = [[Symbol(f'dA{j}{k}') for k in range(3)] for j in range(3)]

    corrections = s3_corrections(A)
    F_tilde_fn  = lambda j, k: make_F_tilde(j, k, A, dA, corrections)
    omega1_fn   = make_omega1(F_tilde_fn, r0, L)
    R_GG        = gamma_gamma_ricci(omega1_fn)
    coeffs      = extract_all(R_GG, A, dA)

Validation::

    from dppu.kk.validator import validate_kk_routes
    validate_kk_routes('t3')    # fast sanity check
    validate_kk_routes('nil3')  # Nil³
    validate_kk_routes('s3')    # S³ (slow)

Modules
-------
field_strength  : F̃, ω⁽¹⁾, topology correction dicts
gamma_gamma     : Γ×Γ shortcut for O(A²) Ricci scalar
full_riemann    : Full Riemann validation route
extractor       : Maxwell / CS / mass² coefficient extraction
validator       : Cross-route validation utility

Author: Muacca
Version: 1.0
"""

from .field_strength import (
    make_F_plain,
    make_F_tilde,
    s3_corrections,
    nil3_corrections,
    t3_corrections,
    make_omega1,
    omega1_to_array,
)

from .gamma_gamma import gamma_gamma_ricci

from .full_riemann import (
    make_omega0,
    full_riemann_scalar_a2,
)

from .extractor import (
    extract_maxwell,
    extract_mass,
    extract_cs,
    extract_all,
)

__all__ = [
    # field_strength
    'make_F_plain',
    'make_F_tilde',
    's3_corrections',
    'nil3_corrections',
    't3_corrections',
    'make_omega1',
    'omega1_to_array',
    # gamma_gamma
    'gamma_gamma_ricci',
    # full_riemann
    'make_omega0',
    'full_riemann_scalar_a2',
    # extractor
    'extract_maxwell',
    'extract_mass',
    'extract_cs',
    'extract_all',
]
