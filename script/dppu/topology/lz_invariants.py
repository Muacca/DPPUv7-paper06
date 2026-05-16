"""LZ-native topology invariants.

This module extracts structure-constant data from the unified topology
engine without relying on paper-specific handoff files.
"""

from __future__ import annotations

from functools import lru_cache
from typing import Any

import sympy as sp

from ..torsion.mode import Mode
from ..torsion.nieh_yan import NyVariant
from .unified import FiberMode, make_engine, normalize_topology_str


A, B, C, U, W = sp.symbols("a b c u v", real=True)
Q = sp.symbols("q", positive=True, real=True)

LZ_PARAMETER_SYMBOLS = (A, B, C, U, W)
CANONICAL_TOPOLOGIES = ("S3", "T3", "Nil3", "Sol3")


def _simplify_expr(expr: Any) -> sp.Expr:
    return sp.factor(sp.cancel(sp.simplify(expr)))


@lru_cache(maxsize=None)
def engine_lz_structure_constants(
    topology: str,
    *,
    signature: str = "lorentzian",
    frame_convention: str = "lz_native",
):
    """Return LZ-native dense structure constants and the engine scale symbol."""
    label = normalize_topology_str(topology)
    engine = make_engine(
        label,
        torsion_mode=Mode.AX,
        ny_variant=NyVariant.FULL,
        signature=signature,
        frame_convention=frame_convention,
        skip_antisymmetry_check=True,
        fiber_mode=FiberMode.NONE,
        enable_squash=False,
        enable_shear=False,
        enable_offdiag_shear=False,
    )
    engine.step_E4_1_setup()
    engine.step_E4_2_metric_and_frame()
    params = engine.data["params"]
    scale = params.get("R", params.get("r"))
    if scale is None:
        raise KeyError(f"scale symbol not found for topology {label}")
    return engine.data["structure_constants"], scale


@lru_cache(maxsize=None)
def extract_lz_five_parameter_specialization(
    topology: str,
    *,
    signature: str = "lorentzian",
    frame_convention: str = "lz_native",
) -> dict[sp.Symbol, sp.Expr]:
    """Extract the five-parameter LZ spatial specialization.

    The returned dictionary maps ``(a,b,c,u,v)`` to the scale-stripped
    components

    ``C^1_23, C^2_31, C^3_12, C^2_12, C^3_13``.
    """
    structure, scale = engine_lz_structure_constants(
        topology,
        signature=signature,
        frame_convention=frame_convention,
    )
    return {
        A: _simplify_expr(structure[1, 2, 3] * scale),
        B: _simplify_expr(structure[2, 3, 1] * scale),
        C: _simplify_expr(structure[3, 1, 2] * scale),
        U: _simplify_expr(structure[2, 1, 2] * scale),
        W: _simplify_expr(structure[3, 1, 3] * scale),
    }


def topology_lz_parameter_table(
    topologies: tuple[str, ...] = CANONICAL_TOPOLOGIES,
    *,
    signature: str = "lorentzian",
    frame_convention: str = "lz_native",
) -> dict[str, dict[sp.Symbol, sp.Expr]]:
    """Return five-parameter LZ specializations for canonical topologies."""
    return {
        normalize_topology_str(topology): extract_lz_five_parameter_specialization(
            topology,
            signature=signature,
            frame_convention=frame_convention,
        )
        for topology in topologies
    }
