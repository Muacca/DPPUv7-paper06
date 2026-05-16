"""
DPPU: Differential Geometry Package for Physics Universe
=========================================================

A modular Python package for Einstein-Cartan gravity with Nieh-Yan term,
designed for publication-quality computations.

Package Structure:
- geometry/    : Metric, volume forms, structure constants
- connection/  : Levi-Civita, contortion, EC connection
- curvature/   : Riemann, Ricci, Hodge dual, Pontryagin
- torsion/     : Mode enum, ansatz, Nieh-Yan density
- action/      : Lagrangian, effective potential, stability, EC action
- topology/    : S3xS1, T3xS1, Nil3xS1 implementations (UnifiedEngine)
- engine/      : Computation pipeline, logging, checkpoints
- utils/       : Shared utilities (epsilon, printing, symbolic, logger)
- kk/          : KK photon effective theory pipeline (v2.1)
                 Two-route computation: Gamma×Gamma shortcut + full Riemann
                 validation.  Supports T³, Nil³, S³ topologies.

Author: Muacca
Version: 2.2
Date: 2026-03
"""

__version__ = "2.2.0"
__author__ = "Muacca"

# Lazy imports for main components
def __getattr__(name):
    """Lazy loading of submodules."""
    if name == "Mode":
        from .torsion.mode import Mode
        return Mode
    elif name == "NyVariant":
        from .torsion.nieh_yan import NyVariant
        return NyVariant
    elif name == "kk":
        from . import kk
        return kk
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
