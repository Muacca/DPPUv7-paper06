"""Spatial Lie-frame curvature helpers."""

from __future__ import annotations

from typing import Any

import sympy as sp

from ..topology.lz_invariants import A, B, C, U, W


def _set_structure_constant(arr, upper: int, lower1: int, lower2: int, value: Any) -> None:
    arr[upper][lower1][lower2] = value
    arr[upper][lower2][lower1] = -value


def spatial_structure_constants_3d(
    params: dict[sp.Symbol, sp.Expr] | None = None,
):
    """Return scale-stripped 3D structure constants for the LZ family."""
    params = params or {A: A, B: B, C: C, U: U, W: W}
    dim = 3
    arr = [[[sp.Integer(0) for _ in range(dim)] for _ in range(dim)] for _ in range(dim)]
    _set_structure_constant(arr, 0, 1, 2, params[A])
    _set_structure_constant(arr, 1, 2, 0, params[B])
    _set_structure_constant(arr, 2, 0, 1, params[C])
    _set_structure_constant(arr, 1, 0, 1, params[U])
    _set_structure_constant(arr, 2, 0, 2, params[W])
    return arr


def scalar_curvature_q2_from_structure_constants(structure_constants) -> sp.Expr:
    """Compute the scale-stripped 3D scalar curvature R^(3) q^2."""
    dim = 3
    gamma = [[[sp.Integer(0) for _ in range(dim)] for _ in range(dim)] for _ in range(dim)]
    for upper in range(dim):
        for i in range(dim):
            for j in range(dim):
                gamma[upper][i][j] = sp.simplify(
                    (
                        structure_constants[upper][i][j]
                        - structure_constants[i][j][upper]
                        + structure_constants[j][upper][i]
                    )
                    / 2
                )

    riemann = [[[[sp.Integer(0) for _ in range(dim)] for _ in range(dim)] for _ in range(dim)] for _ in range(dim)]
    for upper in range(dim):
        for k in range(dim):
            for i in range(dim):
                for j in range(dim):
                    value = sum(
                        gamma[m][j][k] * gamma[upper][i][m]
                        - gamma[m][i][k] * gamma[upper][j][m]
                        for m in range(dim)
                    )
                    value -= sum(
                        structure_constants[m][i][j] * gamma[upper][m][k]
                        for m in range(dim)
                    )
                    riemann[upper][k][i][j] = sp.simplify(value)

    scalar = sum(riemann[i][j][i][j] for i in range(dim) for j in range(dim))
    return sp.factor(sp.cancel(scalar))


def chi_from_lz_parameters(params: dict[sp.Symbol, sp.Expr]) -> sp.Expr:
    """Derive chi = R^(3) q^2 / 6 from LZ five-parameter data."""
    return sp.factor(
        sp.cancel(
            scalar_curvature_q2_from_structure_constants(
                spatial_structure_constants_3d(params)
            )
            / 6
        )
    )
