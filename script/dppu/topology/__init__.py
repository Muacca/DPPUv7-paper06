"""
DPPU topology layer.

Use UnifiedEngine with an explicit DOFConfig.  Historical paper labels and
their configuration equivalents are documented in scripts/README.md.
"""

from .unified import (
    UnifiedEngine,
    DOFConfig,
    TopologyType,
    FiberMode,
    CurvatureSource,
    make_config,
    make_engine,
    normalize_topology,
    normalize_topology_str,
    normalize_curvature_source,
)
from .lz_invariants import (
    A,
    B,
    C,
    U,
    W,
    Q,
    LZ_PARAMETER_SYMBOLS,
    CANONICAL_TOPOLOGIES,
    engine_lz_structure_constants,
    extract_lz_five_parameter_specialization,
    topology_lz_parameter_table,
)
from .base_topology import TopologyEngine
from .s3 import S3Geometry
from .t3 import T3Geometry
from .nil3 import Nil3Geometry
from .sol3 import Sol3Geometry
from .bianchi_ix import (
    bianchi_ix_symbols,
    build_bianchi_ix_lc_weyl,
    build_isotropic_ec_torsion_weyl,
    beta_quadratic_expansion,
)

__all__ = [
    "UnifiedEngine",
    "DOFConfig",
    "TopologyType",
    "FiberMode",
    "CurvatureSource",
    "make_config",
    "make_engine",
    "normalize_topology",
    "normalize_topology_str",
    "normalize_curvature_source",
    "A",
    "B",
    "C",
    "U",
    "W",
    "Q",
    "LZ_PARAMETER_SYMBOLS",
    "CANONICAL_TOPOLOGIES",
    "engine_lz_structure_constants",
    "extract_lz_five_parameter_specialization",
    "topology_lz_parameter_table",
    "TopologyEngine",
    "S3Geometry",
    "T3Geometry",
    "Nil3Geometry",
    "Sol3Geometry",
    "bianchi_ix_symbols",
    "build_bianchi_ix_lc_weyl",
    "build_isotropic_ec_torsion_weyl",
    "beta_quadratic_expansion",
]
