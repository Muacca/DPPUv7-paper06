"""
UnifiedEngine - configurable engine for all geometric degrees of freedom.

The topology layer is configured explicitly with DOFConfig. Historical paper
labels are documented in scripts/README.md; runtime code should use explicit
topology, deformation, and fiber fields.

Usage
-----
    from dppu.topology import make_config, UnifiedEngine, TopologyType, FiberMode
    from dppu.torsion.mode import Mode

    cfg = make_config(
        TopologyType.S3,
        fiber_mode=FiberMode.TWIST,
        isotropic_twist=True,
        torsion_mode=Mode.AX,
    )
    engine = UnifiedEngine(cfg)
    engine.run()
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, Optional, Union

from sympy import symbols, Matrix, S, pi, Rational, sqrt
from sympy.tensor.array import MutableDenseNDimArray

from ..engine.pipeline import BaseFrameEngine
from ..torsion.mode import Mode
from ..torsion.nieh_yan import NyVariant
from ..utils.epsilon import epsilon_3d


# ── Enumerations ──────────────────────────────────────────────────────────────

class TopologyType(Enum):
    """Base topology of the M³×S¹ manifold."""
    S3   = "s3"    # S³×S¹  — SU(2) Lie group × circle
    T3   = "t3"    # T³×S¹  — flat 3-torus × circle
    NIL3 = "nil3"  # Nil³×S¹ — Heisenberg nilmanifold × circle
    SOL3 = "sol3"  # Sol×S¹  — solvable Lie group (Thurston geometry)


class FiberMode(Enum):
    """Deformation type applied to the S¹ fiber direction."""
    NONE    = "none"    # ê³ = L dτ  (unmodified)
    TWIST   = "twist"   # ê³ = L(dτ + ω₁σ⁰ + ω₂σ¹ + ω₃σ²)
    MIXING  = "mixing"  # ê², ê³ rotated by mixing angle δ  (S³ only)
    BOTH    = "both"    # TWIST + MIXING simultaneously (S³ only)
                        # Structure constants computed via full frame rotation
                        # (includes O(δω) cross terms automatically)


class CurvatureSource(Enum):
    """Connection used for curvature-derived diagnostic/action terms."""
    LC = "lc"  # Levi-Civita connection
    EC = "ec"  # Einstein-Cartan/spin connection


# ── DOF Configuration ─────────────────────────────────────────────────────────

@dataclass
class DOFConfig:
    """
    Configuration of active geometric degrees of freedom.

    Parameters
    ----------
    topology : TopologyType
        Base manifold topology (S3, T3, NIL3).
    enable_squash : bool
        Activate axisymmetric squash parameter ε.
        S3: breaks SO(3) → U(1).  Nil3: scales C²₀₁ by (1+ε)^{−4/3}.
    enable_shear : bool
        Activate 0–1 plane shear parameter s (S3 only).
        Further breaks U(1) symmetry of squash (f₀ ≠ f₁ ≠ f₂).
    enable_offdiag_shear : bool
        Activate off-diagonal shear parameters q₃, q₄, q₅ (S3 only).
        Corresponds to spin-2 generators T₃ (01-plane), T₄ (02-plane), T₅ (12-plane).
        Frame: F = F_diag(ε,s) × G(q₃,q₄,q₅).  At q₃=q₄=q₅=0: G=I (backward compat).

        WARNING — Out-of-Memory risk: Running the full pipeline (engine.run())
        with enable_offdiag_shear=True causes OOM even on machines with 64 GB RAM.
        Steps E4.3a and E4.7 (Riemann tensor) produce 5-variable rational products
        of ~3M SymPy ops that exhaust heap memory before completion.
        For off-diagonal Hessian computations, use the single-variable Method 2
        approach (z = exp(q/√2) rational substitution) or a numerical pipeline.
    fiber_mode : FiberMode
        Type of S¹ fiber deformation (NONE | TWIST | MIXING).
    isotropic_twist : bool
        Used only when fiber_mode=TWIST.
        True  → single symbol ω shared across all three twist components.
        False → independent symbols ω₁, ω₂, ω₃.
    torsion_mode : Mode
        AX (η≠0, V=0) | VT (η=0, V≠0) | MX (both).
    ny_variant : NyVariant
        Nieh-Yan term variant: TT | REE | FULL.

    Helper Functions
    ----------------
    make_config(...) -> DOFConfig
        Build this dataclass with explicit topology, deformation, and fiber
        fields.
    """
    topology:           TopologyType = TopologyType.S3
    enable_squash:      bool         = False
    enable_shear:       bool         = False
    enable_offdiag_shear: bool       = False
    fiber_mode:         FiberMode    = FiberMode.NONE
    isotropic_twist:    bool         = True
    torsion_mode:       Mode         = Mode.AX
    ny_variant:         NyVariant    = NyVariant.FULL
    enable_velocity:    bool         = False  # Velocity mode for G metric computation (BOTH mode only)
    skip_antisymmetry_check: bool    = False  # Skip Riemann antisymmetry verification (speedup Option B)
                                              # Enable only when correctness is guaranteed by construction,
                                              # e.g., off-diagonal computations using rationalized variables (z_i).
    weyl_source:        CurvatureSource = CurvatureSource.EC
    riemann_source:     CurvatureSource = CurvatureSource.EC
    signature:          str = "euclidean"       # "euclidean" | "lorentzian"
    frame_convention:   str = "legacy_euclidean" # "legacy_euclidean" | "lz_native"


def normalize_topology(topology: Union[str, TopologyType]) -> TopologyType:
    """Return a TopologyType from enum values or common CLI strings."""
    if isinstance(topology, TopologyType):
        return topology
    if not isinstance(topology, str):
        raise TypeError(f"topology must be str or TopologyType, got {type(topology)}")

    key = topology.strip().lower()
    key = key.replace("_", "").replace("-", "").replace(" ", "")
    key = key.replace("×", "x")
    key = key.replace("xs1", "")
    aliases = {
        "s3": TopologyType.S3,
        "s^3": TopologyType.S3,
        "t3": TopologyType.T3,
        "t^3": TopologyType.T3,
        "nil3": TopologyType.NIL3,
        "nil": TopologyType.NIL3,
        "sol3": TopologyType.SOL3,
        "sol": TopologyType.SOL3,
    }
    try:
        return aliases[key]
    except KeyError as exc:
        raise ValueError(
            f"Unknown topology {topology!r}. Use one of: S3, T3, Nil3, Sol3."
        ) from exc


_TOPOLOGY_LABELS: dict[TopologyType, str] = {
    TopologyType.S3:   "S3",
    TopologyType.T3:   "T3",
    TopologyType.NIL3: "Nil3",
    TopologyType.SOL3: "Sol3",
}


def normalize_topology_str(label: str) -> str:
    """Normalize a topology label to its canonical string form.

    Returns "S3", "T3", "Nil3", or "Sol3".  Accepts the same aliases as
    :func:`normalize_topology`.
    """
    return _TOPOLOGY_LABELS[normalize_topology(label)]


def normalize_curvature_source(
    source: Union[str, CurvatureSource]
) -> CurvatureSource:
    """Return a CurvatureSource from enum values or common CLI strings."""
    if isinstance(source, CurvatureSource):
        return source
    if not isinstance(source, str):
        raise TypeError(
            f"curvature source must be str or CurvatureSource, got {type(source)}"
        )

    key = source.strip().lower()
    aliases = {
        "lc": CurvatureSource.LC,
        "levicivita": CurvatureSource.LC,
        "levi-civita": CurvatureSource.LC,
        "levi_civita": CurvatureSource.LC,
        "ec": CurvatureSource.EC,
        "einsteincartan": CurvatureSource.EC,
        "einstein-cartan": CurvatureSource.EC,
        "spin": CurvatureSource.EC,
    }
    try:
        return aliases[key]
    except KeyError as exc:
        raise ValueError(
            f"Unknown curvature source {source!r}. Use one of: LC, EC."
        ) from exc


def make_config(
    topology: Union[str, TopologyType],
    *,
    enable_squash: bool = False,
    enable_shear: bool = False,
    enable_offdiag_shear: bool = False,
    fiber_mode: FiberMode = FiberMode.NONE,
    isotropic_twist: bool = True,
    torsion_mode: Mode = Mode.AX,
    ny_variant: NyVariant = NyVariant.FULL,
    enable_velocity: bool = False,
    skip_antisymmetry_check: bool = False,
    weyl_source: Union[str, CurvatureSource] = CurvatureSource.EC,
    riemann_source: Union[str, CurvatureSource] = CurvatureSource.EC,
    signature: str = "euclidean",
    frame_convention: str = None,  # None -> auto-derive from signature
) -> DOFConfig:
    """Build a DOFConfig using explicit modern fields.

    Parameters
    ----------
    signature : str
        Metric signature: "euclidean" (default) or "lorentzian".
        Euclidean is the default to preserve existing paper01-04 behaviour.
    frame_convention : str or None
        Frame index convention: "legacy_euclidean" or "lz_native".
        When None (default), derived automatically from signature:
          euclidean  -> "legacy_euclidean"
          lorentzian -> "lz_native"
        An explicit value that conflicts with signature triggers a warning
        but is accepted (recorded in engine.data for auditing).
    """
    import warnings

    _valid_sig = {"euclidean", "lorentzian"}
    _valid_fc  = {"legacy_euclidean", "lz_native"}

    if signature not in _valid_sig:
        raise ValueError(
            f"signature must be one of {_valid_sig}, got {signature!r}"
        )

    # Auto-derive frame_convention when not supplied
    if frame_convention is None:
        frame_convention = "lz_native" if signature == "lorentzian" else "legacy_euclidean"

    if frame_convention not in _valid_fc:
        raise ValueError(
            f"frame_convention must be one of {_valid_fc}, got {frame_convention!r}"
        )

    # Warn on inconsistent combinations
    if signature == "lorentzian" and frame_convention == "legacy_euclidean":
        warnings.warn(
            "signature='lorentzian' with frame_convention='legacy_euclidean' is "
            "inconsistent. Prefer frame_convention='lz_native' for Lorentzian runs.",
            UserWarning,
            stacklevel=2,
        )
    if signature == "euclidean" and frame_convention == "lz_native":
        warnings.warn(
            "signature='euclidean' with frame_convention='lz_native' is inconsistent. "
            "Prefer frame_convention='legacy_euclidean' for Euclidean runs.",
            UserWarning,
            stacklevel=2,
        )

    return DOFConfig(
        topology=normalize_topology(topology),
        enable_squash=enable_squash,
        enable_shear=enable_shear,
        enable_offdiag_shear=enable_offdiag_shear,
        fiber_mode=fiber_mode,
        isotropic_twist=isotropic_twist,
        torsion_mode=torsion_mode,
        ny_variant=ny_variant,
        enable_velocity=enable_velocity,
        skip_antisymmetry_check=skip_antisymmetry_check,
        weyl_source=normalize_curvature_source(weyl_source),
        riemann_source=normalize_curvature_source(riemann_source),
        signature=signature,
        frame_convention=frame_convention,
    )


def make_engine(
    topology: Union[str, TopologyType],
    torsion_mode: Mode = Mode.AX,
    ny_variant: NyVariant = NyVariant.FULL,
    logger=None,
    checkpoint_mgr=None,
    **config_kwargs,
):
    """Build a UnifiedEngine from explicit topology/configuration fields."""
    cfg = make_config(
        topology,
        torsion_mode=torsion_mode,
        ny_variant=ny_variant,
        **config_kwargs,
    )
    return UnifiedEngine(cfg, logger, checkpoint_mgr)


# ── Unified Engine factory ────────────────────────────────────────────────────

class UnifiedEngine:
    """
    Accepts a DOFConfig and returns an instance of the appropriate
    topology implementation class (S3Geometry / T3Geometry / Nil3Geometry).

    Interface::

        engine = UnifiedEngine(cfg, logger)
        engine.run()
        engine.data['potential']     # works as before
        engine.get_riemann_lambdified()
        engine.get_free_params()
        engine.get_config_summary()

    Parameters
    ----------
    config : DOFConfig
        Geometry and torsion configuration (topology / DOFs / mode).
    logger : optional
        Logger instance; NullLogger is used when omitted.
    checkpoint_mgr : optional
        Checkpoint manager for step-level checkpointing.
    """

    def __new__(
        cls,
        config: "DOFConfig",
        logger=None,
        checkpoint_mgr=None,
    ):
        from .s3   import S3Geometry
        from .t3   import T3Geometry
        from .nil3 import Nil3Geometry
        from .sol3 import Sol3Geometry

        if not isinstance(config, DOFConfig):
            raise TypeError(
                f"config must be DOFConfig, got {type(config)}"
            )

        _MAP = {
            TopologyType.S3:   S3Geometry,
            TopologyType.T3:   T3Geometry,
            TopologyType.NIL3: Nil3Geometry,
            TopologyType.SOL3: Sol3Geometry,
        }
        impl_cls = _MAP.get(config.topology)
        if impl_cls is None:
            raise ValueError(
                f"Unknown topology: {config.topology!r}"
            )

        # Delegate to the concrete topology class.
        # Its __init__ (TopologyEngine.__init__) is called normally.
        return impl_cls(config, logger, checkpoint_mgr)

