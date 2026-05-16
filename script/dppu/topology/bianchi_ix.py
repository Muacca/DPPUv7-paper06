"""
Bianchi IX (Squash-Deformed S3) Geometry
==========================================

Provides symbolic computation for the **axisymmetric Bianchi IX**
(squash-deformed S3) geometry in the LZ-native Lorentzian frame::

    e^0 = N dt
    e^1 = r exp( β) σ₁
    e^2 = r exp( β) σ₂
    e^3 = r exp(-2β) σ₃

where ``β(t)`` is the squash parameter and the spatial normalization is
controlled by ``spatial_norm`` (4 for the DPPU S3 convention).

The module exposes two cached builder functions:

* :func:`build_bianchi_ix_lc_weyl` — Levi-Civita Weyl scalar for the
  pure (torsionless) squash geometry.
* :func:`build_isotropic_ec_torsion_weyl` — LC Weyl scalar evaluated at
  β = 0 with EC torsion and time-dependent torsion scalars.

Auxiliary helpers used in Phase2C audits:

* :func:`beta_quadratic_expansion` — small-squash Taylor expansion.
* :func:`bianchi_ix_symbols` — factory for the canonical symbol set.

Symbol conventions
------------------
All symbols are real and declared with SymPy::

    r, rd, rdd   : scale factor and its time derivatives
    N, Nd        : lapse and its time derivative
    beta, bd, bdd: squash parameter and its time derivatives
    alpha        : Weyl-squared coupling
    V3           : spatial volume normalisation (= 2π² for S3)
    kappa        : gravitational coupling
    eta, etad    : axial torsion scalar and its time derivative
    V, Vd        : vector-trace torsion scalar and its time derivative
    theta_NY     : Nieh-Yan coupling

Author: Muacca
"""

from __future__ import annotations

from functools import lru_cache
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple

from sympy import (
    Matrix,
    Rational,
    S,
    cancel,
    diff,
    exp,
    factor,
    series,
    simplify,
    symbols,
)
from sympy.tensor.array import MutableDenseNDimArray

from ..utils.epsilon import epsilon_3d


# ---------------------------------------------------------------------------
# Canonical symbol set
# ---------------------------------------------------------------------------

def bianchi_ix_symbols() -> SimpleNamespace:
    """
    Return the canonical Bianchi IX symbol set as a SimpleNamespace.

    Returns
    -------
    ns : SimpleNamespace
        Attributes: r, rd, rdd, N, Nd, beta, bd, bdd, alpha, V3, kappa,
        eta, etad, V, Vd, theta_NY.
    """
    r, rd, rdd = symbols("r rd rdd", real=True)
    N, Nd = symbols("N Nd", real=True)
    beta, bd, bdd = symbols("beta bd bdd", real=True)
    alpha, V3, kappa = symbols("alpha V3 kappa", real=True)
    eta, etad = symbols("eta etad", real=True)
    V, Vd = symbols("V Vd", real=True)
    theta_NY = symbols("theta_NY", real=True)
    return SimpleNamespace(
        r=r, rd=rd, rdd=rdd,
        N=N, Nd=Nd,
        beta=beta, bd=bd, bdd=bdd,
        alpha=alpha, V3=V3, kappa=kappa,
        eta=eta, etad=etad,
        V=V, Vd=Vd,
        theta_NY=theta_NY,
    )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

_DIM = 4
_METRIC = Matrix.diag(-1, 1, 1, 1)
_METRIC_INV = Matrix.diag(-1, 1, 1, 1)


def _lower_c(C: MutableDenseNDimArray, a: int, b: int, c: int) -> Any:
    return sum(_METRIC[a, d] * C[d, b, c] for d in range(_DIM))


def _structure_constants(sy: SimpleNamespace, spatial_norm: int = 4) -> MutableDenseNDimArray:
    """Bianchi IX squash-deformed S3 structure constants."""
    H = [
        (sy.rd / sy.r + sy.bd) / sy.N,
        (sy.rd / sy.r + sy.bd) / sy.N,
        (sy.rd / sy.r - 2 * sy.bd) / sy.N,
    ]
    lam = [
        spatial_norm * exp(2 * sy.beta) / sy.r,
        spatial_norm * exp(2 * sy.beta) / sy.r,
        spatial_norm * exp(-4 * sy.beta) / sy.r,
    ]
    C = MutableDenseNDimArray.zeros(_DIM, _DIM, _DIM)
    for i in range(3):
        C[i + 1, 0, i + 1] = H[i]
        C[i + 1, i + 1, 0] = -H[i]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                eps = epsilon_3d(i, j, k)
                if eps:
                    C[i + 1, j + 1, k + 1] = lam[i] * eps
    return C


def _koszul_connection(C: MutableDenseNDimArray) -> MutableDenseNDimArray:
    """Metric-aware Levi-Civita connection for Lorentzian signature."""
    Gamma = MutableDenseNDimArray.zeros(_DIM, _DIM, _DIM)
    for a in range(_DIM):
        for b in range(_DIM):
            for c in range(_DIM):
                val = S.Zero
                for e in range(_DIM):
                    gamma_low = (
                        _lower_c(C, c, b, e)
                        - _lower_c(C, b, e, c)
                        - _lower_c(C, e, c, b)
                    ) / 2
                    val += _METRIC_INV[a, e] * gamma_low
                Gamma[a, b, c] = cancel(val)
    return Gamma


def _frame_d0(
    expr,
    sy: SimpleNamespace,
    extra: Optional[List[Tuple]] = None,
) -> Any:
    """
    Frame time-derivative :math:`E_0(f) = N^{-1} \\dot f`.

    Handles ``r, rd, rdd, N, Nd, beta, bd, bdd`` plus any additional
    ``(symbol, time_derivative)`` pairs supplied via *extra*.
    """
    pairs = [
        (sy.r, sy.rd),
        (sy.rd, sy.rdd),
        (sy.N, sy.Nd),
        (sy.beta, sy.bd),
        (sy.bd, sy.bdd),
    ]
    if extra:
        pairs.extend(extra)
    val = sum(diff(expr, q) * qd for q, qd in pairs)
    return cancel(val / sy.N)


def _frame_deriv(c: int, expr, sy: SimpleNamespace, extra=None) -> Any:
    return _frame_d0(expr, sy, extra=extra) if c == 0 else S.Zero


def _riemann_from_connection(
    Gamma: MutableDenseNDimArray,
    C: MutableDenseNDimArray,
    sy: SimpleNamespace,
    extra_derivatives=None,
) -> MutableDenseNDimArray:
    """Riemann tensor :math:`R^a{}_{bcd}` from connection coefficients."""
    Riemann = MutableDenseNDimArray.zeros(_DIM, _DIM, _DIM, _DIM)
    for a in range(_DIM):
        for b in range(_DIM):
            for c in range(_DIM):
                for d in range(_DIM):
                    term = _frame_deriv(
                        c, Gamma[a, b, d], sy, extra=extra_derivatives
                    ) - _frame_deriv(
                        d, Gamma[a, b, c], sy, extra=extra_derivatives
                    )
                    for e in range(_DIM):
                        term += Gamma[e, b, d] * Gamma[a, e, c]
                        term -= Gamma[e, b, c] * Gamma[a, e, d]
                        term += C[e, c, d] * Gamma[a, b, e]
                    Riemann[a, b, c, d] = cancel(term)
    return Riemann


def _lower_first_index(Riemann: MutableDenseNDimArray) -> MutableDenseNDimArray:
    R_abcd = MutableDenseNDimArray.zeros(_DIM, _DIM, _DIM, _DIM)
    for a in range(_DIM):
        for b in range(_DIM):
            for c in range(_DIM):
                for d in range(_DIM):
                    R_abcd[a, b, c, d] = sum(
                        _METRIC[a, e] * Riemann[e, b, c, d] for e in range(_DIM)
                    )
    return R_abcd


def _ricci_tensor(Riemann: MutableDenseNDimArray) -> Matrix:
    Ricci = Matrix.zeros(_DIM, _DIM)
    for b in range(_DIM):
        for d in range(_DIM):
            Ricci[b, d] = cancel(sum(Riemann[a, b, a, d] for a in range(_DIM)))
    return Ricci


def _ricci_scalar(Ricci: Matrix) -> Any:
    return cancel(
        sum(_METRIC_INV[a, b] * Ricci[a, b] for a in range(_DIM) for b in range(_DIM))
    )


def _weyl_tensor(
    R_abcd: MutableDenseNDimArray,
    Ricci: Matrix,
    R_scalar,
) -> MutableDenseNDimArray:
    C_abcd = MutableDenseNDimArray.zeros(_DIM, _DIM, _DIM, _DIM)
    for a in range(_DIM):
        for b in range(_DIM):
            for c in range(_DIM):
                for d in range(_DIM):
                    term = R_abcd[a, b, c, d]
                    term -= Rational(1, 2) * (
                        _METRIC[a, c] * Ricci[b, d]
                        - _METRIC[a, d] * Ricci[b, c]
                        - _METRIC[b, c] * Ricci[a, d]
                        + _METRIC[b, d] * Ricci[a, c]
                    )
                    term += Rational(1, 6) * R_scalar * (
                        _METRIC[a, c] * _METRIC[b, d]
                        - _METRIC[a, d] * _METRIC[b, c]
                    )
                    C_abcd[a, b, c, d] = cancel(term)
    return C_abcd


def _weyl_scalar(C_abcd: MutableDenseNDimArray) -> Any:
    C2 = S.Zero
    for a in range(_DIM):
        for b in range(_DIM):
            for c in range(_DIM):
                for d in range(_DIM):
                    C2 += (
                        C_abcd[a, b, c, d]
                        * _METRIC_INV[a, a]
                        * _METRIC_INV[b, b]
                        * _METRIC_INV[c, c]
                        * _METRIC_INV[d, d]
                        * C_abcd[a, b, c, d]
                    )
    return cancel(C2)


def _lz_native_torsion(mode: str, sy: SimpleNamespace) -> MutableDenseNDimArray:
    """Build LZ-native torsion tensor :math:`T_{abc}` at isotropic limit."""
    T = MutableDenseNDimArray.zeros(_DIM, _DIM, _DIM)
    if mode in {"AX", "MX"}:
        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    eps = epsilon_3d(i - 1, j - 1, k - 1)
                    if eps:
                        T[i, j, k] = T[i, j, k] + 2 * sy.eta * eps / sy.r
    if mode in {"VT", "MX"}:
        trace = [sy.V, S.Zero, S.Zero, S.Zero]
        for a in range(_DIM):
            for b in range(_DIM):
                for c in range(_DIM):
                    val = (_METRIC[a, c] * trace[b] - _METRIC[a, b] * trace[c]) / 3
                    if val != S.Zero:
                        T[a, b, c] = T[a, b, c] + val
    return T


def _contortion_from_torsion(T: MutableDenseNDimArray) -> MutableDenseNDimArray:
    K_low = MutableDenseNDimArray.zeros(_DIM, _DIM, _DIM)
    K = MutableDenseNDimArray.zeros(_DIM, _DIM, _DIM)
    for a in range(_DIM):
        for b in range(_DIM):
            for c in range(_DIM):
                K_low[a, b, c] = cancel((T[a, b, c] + T[b, c, a] - T[c, a, b]) / 2)
    for a in range(_DIM):
        for b in range(_DIM):
            for c in range(_DIM):
                K[a, b, c] = cancel(
                    sum(_METRIC_INV[a, d] * K_low[d, b, c] for d in range(_DIM))
                )
    return K


# ---------------------------------------------------------------------------
# Public cached builders
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def build_bianchi_ix_lc_weyl(spatial_norm: int = 4) -> Dict:
    """
    Build the LC Weyl scalar for the axisymmetric Bianchi IX geometry.

    The result is cached by *spatial_norm*.

    Returns
    -------
    dict with keys:
      ``symbols``, ``structure_constants``, ``connection``,
      ``riemann``, ``riemann_abcd``, ``ricci``, ``ricci_scalar``,
      ``weyl_tensor``, ``weyl_scalar``, ``spatial_norm``, ``chi``.
    """
    sy = bianchi_ix_symbols()
    C = _structure_constants(sy, spatial_norm=spatial_norm)
    Gamma = _koszul_connection(C)
    Riemann = _riemann_from_connection(Gamma, C, sy)
    R_abcd = _lower_first_index(Riemann)
    Ricci = _ricci_tensor(Riemann)
    R_scalar = _ricci_scalar(Ricci)
    C_abcd = _weyl_tensor(R_abcd, Ricci, R_scalar)
    C2 = _weyl_scalar(C_abcd)
    return {
        "symbols": sy,
        "structure_constants": C,
        "connection": Gamma,
        "riemann": Riemann,
        "riemann_abcd": R_abcd,
        "ricci": Ricci,
        "ricci_scalar": R_scalar,
        "weyl_tensor": C_abcd,
        "weyl_scalar": C2,
        "spatial_norm": spatial_norm,
        "chi": Rational(spatial_norm, 2) ** 2,
    }


@lru_cache(maxsize=None)
def build_isotropic_ec_torsion_weyl(mode: str, spatial_norm: int = 4) -> Dict:
    """
    Build the EC Weyl scalar in the isotropic limit (β = 0) with
    time-dependent torsion scalars η(t), V(t).

    The isotropic limit is imposed by substituting β = β̇ = β̈ = 0
    in the structure constants before adding EC torsion.

    Parameters
    ----------
    mode:
        Torsion mode — one of ``"AX"``, ``"VT"``, ``"MX"``.
    spatial_norm:
        S3 normalization factor (default 4).

    Returns
    -------
    dict with keys:
      ``symbols``, ``contortion``, ``riemann``, ``riemann_abcd``,
      ``ricci``, ``ricci_scalar``, ``weyl_tensor``, ``weyl_scalar``,
      ``mode``.
    """
    mode = mode.upper()
    if mode not in {"AX", "VT", "MX"}:
        raise ValueError(f"mode must be AX, VT, or MX; got {mode!r}")

    sy = bianchi_ix_symbols()
    C_full = _structure_constants(sy, spatial_norm=spatial_norm)
    iso_subs = {sy.beta: S.Zero, sy.bd: S.Zero, sy.bdd: S.Zero}
    C = MutableDenseNDimArray.zeros(_DIM, _DIM, _DIM)
    for a in range(_DIM):
        for b in range(_DIM):
            for c in range(_DIM):
                C[a, b, c] = simplify(S(C_full[a, b, c]).subs(iso_subs))

    Gamma_lc = _koszul_connection(C)
    K = _contortion_from_torsion(_lz_native_torsion(mode, sy))
    Gamma_ec = MutableDenseNDimArray.zeros(_DIM, _DIM, _DIM)
    for a in range(_DIM):
        for b in range(_DIM):
            for c in range(_DIM):
                Gamma_ec[a, b, c] = cancel(Gamma_lc[a, b, c] + K[a, b, c])

    extra = [(sy.eta, sy.etad), (sy.V, sy.Vd)]
    Riemann = _riemann_from_connection(Gamma_ec, C, sy, extra_derivatives=extra)
    R_abcd = _lower_first_index(Riemann)
    Ricci = _ricci_tensor(Riemann)
    R_scalar = _ricci_scalar(Ricci)
    C_abcd = _weyl_tensor(R_abcd, Ricci, R_scalar)
    C2 = factor(_weyl_scalar(C_abcd))
    return {
        "symbols": sy,
        "contortion": K,
        "riemann": Riemann,
        "riemann_abcd": R_abcd,
        "ricci": Ricci,
        "ricci_scalar": R_scalar,
        "weyl_tensor": C_abcd,
        "weyl_scalar": C2,
        "mode": mode,
    }


# ---------------------------------------------------------------------------
# Squash expansion helper
# ---------------------------------------------------------------------------

def beta_quadratic_expansion(expr) -> Dict:
    """
    Expand *expr* to O(ε²) with β → ε x, β̇ → ε ẋ, β̈ → ε ẍ.

    Parameters
    ----------
    expr:
        SymPy expression depending on the standard Bianchi IX symbols.

    Returns
    -------
    dict with keys ``"c0"``, ``"c1"``, ``"c2"`` for the constant,
    linear, and quadratic coefficients in ε.
    """
    sy = bianchi_ix_symbols()
    eps, x, xd, xdd = symbols("eps x xd xdd", real=True)
    scaled = expr.subs(
        {sy.beta: eps * x, sy.bd: eps * xd, sy.bdd: eps * xdd}
    )
    poly = series(scaled, eps, 0, 3).removeO().expand()
    back = {x: sy.beta, xd: sy.bd, xdd: sy.bdd}
    return {
        "c0": simplify(poly.coeff(eps, 0).subs(back)),
        "c1": simplify(poly.coeff(eps, 1).subs(back)),
        "c2": factor(simplify(poly.coeff(eps, 2).subs(back))),
    }
