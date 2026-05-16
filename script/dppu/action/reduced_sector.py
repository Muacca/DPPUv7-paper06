"""Reduced-sector extraction helpers for EC+NY topology engines."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any

import sympy as sp

from ..curvature.spatial_lie import chi_from_lz_parameters
from ..topology.lz_invariants import (
    CANONICAL_TOPOLOGIES,
    Q,
    engine_lz_structure_constants,
    extract_lz_five_parameter_specialization,
)
from ..topology.unified import FiberMode, make_engine, normalize_topology_str
from ..torsion.mode import Mode, normalize_mode
from ..torsion.nieh_yan import NyVariant


QDOT = sp.symbols("qdot", real=True)
N = sp.symbols("N", positive=True, real=True)
ETA = sp.symbols("eta", real=True)
V_TORSION = sp.symbols("V", real=True)
KAPPA = sp.symbols("kappa", positive=True, real=True)
THETA = sp.symbols("theta_NY", real=True)

TOPOLOGIES = CANONICAL_TOPOLOGIES
MODES = ("EH", "AX", "VT", "MX")

_MODE_ENUM = {
    "AX": Mode.AX,
    "VT": Mode.VT,
    "MX": Mode.MX,
}


@dataclass(frozen=True)
class StaticBranchResult:
    """Engine-derived static reduced branch data."""

    topology: str
    mode: str
    F_branch: sp.Expr
    chi: sp.Expr
    scale_symbol: sp.Symbol
    potential: sp.Expr | None
    substitutions: dict[Any, Any]


def _simplify_expr(expr: Any) -> sp.Expr:
    return sp.factor(sp.cancel(sp.simplify(expr)))


def normalize_reduced_mode(mode: str | Mode) -> str:
    """Normalize ``EH`` or torsion mode labels to canonical strings."""
    value = getattr(mode, "value", mode)
    key = str(value).strip().upper()
    if key == "EH":
        return "EH"
    return normalize_mode(key)


@lru_cache(maxsize=None)
def chi_for_topology(topology: str) -> sp.Expr:
    """Return chi for a topology by LZ structure constants and 3D curvature."""
    return chi_from_lz_parameters(extract_lz_five_parameter_specialization(topology))


def _make_lz_engine(topology: str, mode: str):
    return make_engine(
        normalize_topology_str(topology),
        torsion_mode=_MODE_ENUM[mode],
        ny_variant=NyVariant.FULL,
        signature="lorentzian",
        frame_convention="lz_native",
        skip_antisymmetry_check=True,
        fiber_mode=FiberMode.NONE,
        enable_squash=False,
        enable_shear=False,
        enable_offdiag_shear=False,
    )


@lru_cache(maxsize=None)
def derive_static_branch_result(topology: str, mode: str | Mode) -> StaticBranchResult:
    """Derive the isotropic static reduced branch function from the engine."""
    label = normalize_topology_str(topology)
    mode_label = normalize_reduced_mode(mode)
    chi_value = chi_for_topology(label)

    if mode_label == "EH":
        _, scale_symbol = engine_lz_structure_constants(label)
        return StaticBranchResult(
            topology=label,
            mode=mode_label,
            F_branch=chi_value,
            chi=chi_value,
            scale_symbol=scale_symbol,
            potential=None,
            substitutions={scale_symbol: Q},
        )

    engine = _make_lz_engine(label, mode_label)
    engine.run()
    params = engine.data["params"]
    potential = engine.data["potential"].subs(params["alpha"], sp.S.Zero)
    scale_engine = params.get("R", params.get("r"))
    if scale_engine is None:
        raise KeyError(f"scale symbol not found for topology {label}")

    static_base_volume = sp.simplify(engine.data["total_volume"] / (params["L"] * scale_engine**3))
    action_prefactor = sp.simplify(3 * static_base_volume / params["kappa"] ** 2)
    f_engine = sp.simplify((-potential / params["L"]) / (action_prefactor * scale_engine))

    substitutions = {
        params["kappa"]: KAPPA,
        params["theta_NY"]: THETA,
        scale_engine: Q,
    }
    if params.get("r") is not None:
        substitutions[params["r"]] = Q
    if params.get("R") is not None:
        substitutions[params["R"]] = Q
    if params["eta"] != sp.S.Zero:
        substitutions[params["eta"]] = ETA
    if params["V"] != sp.S.Zero:
        substitutions[params["V"]] = V_TORSION

    return StaticBranchResult(
        topology=label,
        mode=mode_label,
        F_branch=_simplify_expr(f_engine.subs(substitutions)),
        chi=chi_value,
        scale_symbol=scale_engine,
        potential=potential,
        substitutions=substitutions,
    )


def derive_static_branch_function(topology: str, mode: str | Mode) -> sp.Expr:
    """Return only the engine-derived static reduced branch function F."""
    return derive_static_branch_result(topology, mode).F_branch


def auxiliary_variables(expr: sp.Expr, candidates: tuple[sp.Symbol, ...] = (ETA, V_TORSION)) -> list[sp.Symbol]:
    """Return auxiliary variables from ``candidates`` that occur in ``expr``."""
    return [variable for variable in candidates if expr.has(variable)]


def solve_auxiliary_shell(
    expr: sp.Expr,
    aux_vars: list[sp.Symbol] | tuple[sp.Symbol, ...] | None = None,
) -> dict[str, object]:
    """Solve auxiliary equations and compute Hessian data."""
    variables = list(aux_vars) if aux_vars is not None else auxiliary_variables(expr)
    if not variables:
        return {
            "aux_vars": [],
            "equations": [],
            "solution": {},
            "solution_status": "NO_AUXILIARY_VARIABLES",
            "hessian": None,
            "hessian_det": None,
            "hessian_status": "NOT_APPLICABLE",
            "on_shell": sp.factor(sp.cancel(expr)),
        }

    equations = [sp.diff(expr, variable) for variable in variables]
    solutions = sp.solve(equations, variables, dict=True)
    solution = solutions[0] if len(solutions) == 1 else {}
    hessian = sp.Matrix([[sp.diff(expr, x, y) for y in variables] for x in variables])
    hessian_det = sp.factor(sp.cancel(hessian.det()))
    on_shell = _simplify_expr(expr.subs(solution)) if solution else sp.nan
    if len(solutions) != 1:
        solution_status = "NON_UNIQUE_OR_UNSOLVED"
    elif hessian_det == 0:
        solution_status = "DEGENERATE_HESSIAN"
    else:
        solution_status = "UNIQUE_NONDEGENERATE"

    return {
        "aux_vars": variables,
        "equations": equations,
        "solution": solution,
        "solution_status": solution_status,
        "hessian": hessian,
        "hessian_det": hessian_det,
        "hessian_status": "NONZERO_UNDER_Q_POSITIVE_DOMAIN" if hessian_det != 0 else "ZERO",
        "on_shell": on_shell,
    }
