"""
DPPU Curvature Layer
====================

Provides curvature-related components for Einstein-Cartan gravity.

Modules:
- riemann: Riemann tensor computation and antisymmetry verification
- ricci: Ricci tensor and scalar
- hodge: Hodge dual operator
- pontryagin: Pontryagin density P = <R, *R>
- sd_extension: SDExtensionMixin for adding SD methods to engines
- sd_diagnostics: CurvatureSDDiagnostics for SD/ASD analysis
- spatial_lie: 3D Lie-frame scalar curvature helpers
"""

from .riemann import RiemannAntisymmetryError, verify_antisymmetry_strict
from .ricci import compute_ricci_scalar, compute_ricci_tensor
from .hodge import compute_hodge_dual
from .pontryagin import compute_P_from_riemann, get_riemann_numeric
from .sd_extension import SDExtensionMixin
from .sd_diagnostics import CurvatureSDDiagnostics

__all__ = [
    'RiemannAntisymmetryError',
    'verify_antisymmetry_strict',
    'compute_ricci_scalar',
    'compute_ricci_tensor',
    'compute_hodge_dual',
    'compute_P_from_riemann',
    'get_riemann_numeric',
    'SDExtensionMixin',
    'CurvatureSDDiagnostics',
]
