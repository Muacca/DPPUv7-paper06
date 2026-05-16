"""
Nieh-Yan Form Engine (LZ-Native Lorentzian S3)
================================================

Provides :class:`NiehYanFormEngine` — a self-contained suite for
computing Nieh-Yan boundary forms and related diagnostics using
exterior form algebra on the LZ-native Lorentzian S3 coframe::

    e^0 = N dt
    e^i = r σ^i   (i = 1,2,3)

with structure constants :math:`C^i_{jk} = (4/r)\,\\varepsilon_{ijk}`
(DPPU S3 convention with ``spatial_norm=4``).

Usage::

    from dppu.forms.nieh_yan import NiehYanFormEngine

    engine = NiehYanFormEngine(symbols, spatial_norm=4)
    audit  = engine.ny_form_audit("AX")

The *symbols* argument must expose the following attributes::

    r, rd, rdd, N, Nd, eta, etad, V, Vd, V3, kappa, theta_NY

(attribute access, e.g. a dataclass or SimpleNamespace).

Author: Muacca
"""

from __future__ import annotations

from typing import Dict, Iterable, Tuple

from sympy import S, Matrix, cancel, simplify
from sympy.tensor.array import MutableDenseNDimArray

from .exterior_algebra import (
    Form,
    add_forms,
    basis,
    clean_form,
    coefficient,
    form_to_str,
    scale_form,
    wedge,
)
from ..utils.epsilon import epsilon_3d


class NiehYanFormEngine:
    """
    Nieh-Yan form engine for LZ-native Lorentzian S3.

    All physics functions are methods of this class so that the symbol
    set is stored once and reused consistently.

    Parameters
    ----------
    symbols:
        An object with attributes ``r``, ``rd``, ``rdd``, ``N``, ``Nd``,
        ``eta``, ``etad``, ``V``, ``Vd``, ``V3``, ``kappa``, ``theta_NY``.
    spatial_norm:
        S3 normalization factor (default 4, DPPU convention).
    """

    DIM: int = 4
    METRIC: Matrix = Matrix.diag(-1, 1, 1, 1)
    METRIC_INV: Matrix = Matrix.diag(-1, 1, 1, 1)

    def __init__(self, symbols, spatial_norm: int = 4) -> None:
        self.sy = symbols
        self.spatial_norm = spatial_norm

    # ------------------------------------------------------------------
    # Structure constants
    # ------------------------------------------------------------------

    def structure_constants(self) -> MutableDenseNDimArray:
        """
        Build isotropic LZ-native S3 structure constants.

        Returns
        -------
        C : MutableDenseNDimArray
            :math:`C^a_{bc}` with
            :math:`C^i_{0i} = \\dot r / (N r)` and
            :math:`C^i_{jk} = (4/r)\\,\\varepsilon_{ijk}`.
        """
        sy = self.sy
        C = MutableDenseNDimArray.zeros(self.DIM, self.DIM, self.DIM)
        H = sy.rd / (sy.N * sy.r)
        for i in range(1, 4):
            C[i, 0, i] = H
            C[i, i, 0] = -H
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    eps = epsilon_3d(i, j, k)
                    if eps:
                        C[i + 1, j + 1, k + 1] = self.spatial_norm * eps / sy.r
        return C

    # ------------------------------------------------------------------
    # Exterior derivative infrastructure
    # ------------------------------------------------------------------

    def _de_basis_forms(self, C: MutableDenseNDimArray) -> Tuple[Form, ...]:
        """
        Return ``de^a = (1/2) C^a_{bc} e^b ∧ e^c`` for each frame index *a*.
        """
        forms = []
        for a in range(self.DIM):
            item: Form = {}
            for b in range(self.DIM):
                for c in range(self.DIM):
                    val = C[a, b, c] / 2
                    if val != S.Zero:
                        item = add_forms(
                            item, scale_form(val, wedge(basis(b), basis(c)))
                        )
            forms.append(item)
        return tuple(forms)

    def frame_time_derivative(self, expr) -> object:
        """
        Frame time-derivative :math:`E_0(f) = N^{-1} \\dot f`.

        Accounts for ``r, rd, rdd, N, Nd, eta, etad, V, Vd``.
        """
        sy = self.sy
        total = (
            expr.diff(sy.r) * sy.rd
            + expr.diff(sy.rd) * sy.rdd
            + expr.diff(sy.N) * sy.Nd
            + expr.diff(sy.eta) * sy.etad
            + expr.diff(sy.V) * sy.Vd
        )
        return simplify(total / sy.N)

    def exterior_derivative(
        self, form: Form, C: MutableDenseNDimArray | None = None
    ) -> Form:
        """
        Compute the exterior derivative ``d(form)`` using the Cartan
        structure equation :math:`de^a = (1/2) C^a_{bc} e^b \\wedge e^c`.
        """
        if C is None:
            C = self.structure_constants()
        de = self._de_basis_forms(C)
        out: Form = {}
        for idx, coeff in form.items():
            if coeff != S.Zero:
                out = add_forms(
                    out,
                    wedge(
                        {(0,): self.frame_time_derivative(coeff)},
                        {idx: S.One},
                    ),
                )
            for pos, a in enumerate(idx):
                prefix: Form = {idx[:pos]: S.One}
                suffix: Form = {idx[pos + 1:]: S.One}
                term = wedge(wedge(prefix, de[a]), suffix)
                if pos % 2:
                    term = scale_form(-S.One, term)
                out = add_forms(out, scale_form(coeff, term))
        return clean_form(out)

    # ------------------------------------------------------------------
    # Torsion
    # ------------------------------------------------------------------

    def torsion_tensor_low(self, mode: str) -> MutableDenseNDimArray:
        """
        Build the lowered torsion tensor :math:`T_{abc}` for *mode*.

        ``mode`` must be one of ``"AX"``, ``"VT"``, ``"MX"``.
        """
        sy = self.sy
        mode = mode.upper()
        T = MutableDenseNDimArray.zeros(self.DIM, self.DIM, self.DIM)

        if mode in {"AX", "MX"}:
            for i in range(1, 4):
                for j in range(1, 4):
                    for k in range(1, 4):
                        eps = epsilon_3d(i - 1, j - 1, k - 1)
                        if eps:
                            T[i, j, k] = T[i, j, k] + 2 * sy.eta * eps / sy.r

        if mode in {"VT", "MX"}:
            trace = [sy.V, S.Zero, S.Zero, S.Zero]
            for a in range(self.DIM):
                for b in range(self.DIM):
                    for c in range(self.DIM):
                        val = (
                            self.METRIC[a, c] * trace[b]
                            - self.METRIC[a, b] * trace[c]
                        ) / 3
                        if val != S.Zero:
                            T[a, b, c] = T[a, b, c] + val
        return T

    def torsion_forms_low(self, mode: str) -> Tuple[Form, ...]:
        """Return lowered torsion 2-forms :math:`T_a = (1/2) T_{abc} e^b \\wedge e^c`."""
        T = self.torsion_tensor_low(mode)
        forms = []
        for a in range(self.DIM):
            form: Form = {}
            for b in range(self.DIM):
                for c in range(self.DIM):
                    val = T[a, b, c] / 2
                    if val != S.Zero:
                        form = add_forms(
                            form, scale_form(val, wedge(basis(b), basis(c)))
                        )
            forms.append(form)
        return tuple(forms)

    def torsion_forms_up(self, mode: str) -> Tuple[Form, ...]:
        """Return raised torsion 2-forms :math:`T^a = (1/2) T^a_{bc} e^b \\wedge e^c`."""
        T = self.torsion_tensor_low(mode)
        forms = []
        for a in range(self.DIM):
            form: Form = {}
            for b in range(self.DIM):
                for c in range(self.DIM):
                    val = (
                        sum(self.METRIC_INV[a, d] * T[d, b, c] for d in range(self.DIM))
                        / 2
                    )
                    if val != S.Zero:
                        form = add_forms(
                            form, scale_form(val, wedge(basis(b), basis(c)))
                        )
            forms.append(form)
        return tuple(forms)

    # ------------------------------------------------------------------
    # Nieh-Yan 3-form Q and its exterior derivative
    # ------------------------------------------------------------------

    def q_nieh_yan(self, mode: str) -> Form:
        """Return the Nieh-Yan 3-form :math:`Q = e^a \\wedge T_a`."""
        T_low = self.torsion_forms_low(mode)
        qform: Form = {}
        for a in range(self.DIM):
            qform = add_forms(qform, wedge(basis(a), T_low[a]))
        return qform

    def integrated_q_boundary(self, mode: str) -> object:
        """
        Return the integrated boundary functional
        :math:`V_3 \\cdot q_{123} \\cdot r^3`.
        """
        sy = self.sy
        qcoeff = coefficient(self.q_nieh_yan(mode), (1, 2, 3))
        return simplify(sy.V3 * qcoeff * sy.r**3)

    def dq_reduced_density(self, mode: str) -> object:
        """
        Return the reduced bulk density from :math:`dQ`.

        Multiplied by :math:`V_3 N r^3` to give the Lagrangian density.
        """
        sy = self.sy
        dq_coeff = coefficient(self.exterior_derivative(self.q_nieh_yan(mode)))
        return simplify(sy.V3 * sy.N * sy.r**3 * dq_coeff)

    # ------------------------------------------------------------------
    # Auxiliary forms for NY identity
    # ------------------------------------------------------------------

    def tt_form(self, mode: str) -> Form:
        """Return :math:`T^a \\wedge T_a`."""
        up = self.torsion_forms_up(mode)
        low = self.torsion_forms_low(mode)
        out: Form = {}
        for a in range(self.DIM):
            out = add_forms(out, wedge(up[a], low[a]))
        return out

    # ------------------------------------------------------------------
    # EC connection and curvature 2-forms
    # ------------------------------------------------------------------

    def _contortion(self, T: MutableDenseNDimArray) -> MutableDenseNDimArray:
        """Compute contortion :math:`K^a_{bc}` from all-lowered torsion."""
        K_low = MutableDenseNDimArray.zeros(self.DIM, self.DIM, self.DIM)
        K = MutableDenseNDimArray.zeros(self.DIM, self.DIM, self.DIM)
        for a in range(self.DIM):
            for b in range(self.DIM):
                for c in range(self.DIM):
                    K_low[a, b, c] = cancel(
                        (T[a, b, c] + T[b, c, a] - T[c, a, b]) / 2
                    )
        for a in range(self.DIM):
            for b in range(self.DIM):
                for c in range(self.DIM):
                    K[a, b, c] = cancel(
                        sum(
                            self.METRIC_INV[a, d] * K_low[d, b, c]
                            for d in range(self.DIM)
                        )
                    )
        return K

    def _koszul_connection(
        self, C: MutableDenseNDimArray
    ) -> MutableDenseNDimArray:
        """
        Compute the Levi-Civita connection via the metric-aware Koszul formula.

        This version is valid for Lorentzian signature (metric ≠ identity).
        """
        Gamma = MutableDenseNDimArray.zeros(self.DIM, self.DIM, self.DIM)
        for a in range(self.DIM):
            for b in range(self.DIM):
                for c in range(self.DIM):
                    val = S.Zero
                    for e in range(self.DIM):
                        gamma_low = (
                            sum(
                                self.METRIC[c, d] * C[d, b, e]
                                for d in range(self.DIM)
                            )
                            - sum(
                                self.METRIC[b, d] * C[d, e, c]
                                for d in range(self.DIM)
                            )
                            - sum(
                                self.METRIC[e, d] * C[d, c, b]
                                for d in range(self.DIM)
                            )
                        ) / 2
                        val += self.METRIC_INV[a, e] * gamma_low
                    Gamma[a, b, c] = cancel(val)
        return Gamma

    def _riemann_from_connection(
        self,
        Gamma: MutableDenseNDimArray,
        C: MutableDenseNDimArray,
    ) -> MutableDenseNDimArray:
        """Compute :math:`R^a{}_{bcd}` from connection coefficients."""
        Riemann = MutableDenseNDimArray.zeros(
            self.DIM, self.DIM, self.DIM, self.DIM
        )
        for a in range(self.DIM):
            for b in range(self.DIM):
                for c in range(self.DIM):
                    for d in range(self.DIM):
                        term = (
                            self._frame_d0_conn(c, Gamma[a, b, d])
                            - self._frame_d0_conn(d, Gamma[a, b, c])
                        )
                        for e in range(self.DIM):
                            term += Gamma[e, b, d] * Gamma[a, e, c]
                            term -= Gamma[e, b, c] * Gamma[a, e, d]
                            term += C[e, c, d] * Gamma[a, b, e]
                        Riemann[a, b, c, d] = cancel(term)
        return Riemann

    def _frame_d0_conn(self, c: int, expr) -> object:
        """Apply :math:`E_c` (non-zero only for c=0 in cosmological ansatz)."""
        return self.frame_time_derivative(expr) if c == 0 else S.Zero

    def _lower_first_index(
        self, Riemann: MutableDenseNDimArray
    ) -> MutableDenseNDimArray:
        """Lower first index: :math:`R_{abcd} = \\eta_{ae} R^e{}_{bcd}`."""
        R_abcd = MutableDenseNDimArray.zeros(
            self.DIM, self.DIM, self.DIM, self.DIM
        )
        for a in range(self.DIM):
            for b in range(self.DIM):
                for c in range(self.DIM):
                    for d in range(self.DIM):
                        R_abcd[a, b, c, d] = sum(
                            self.METRIC[a, e] * Riemann[e, b, c, d]
                            for e in range(self.DIM)
                        )
        return R_abcd

    def curvature_forms_low(self, mode: str) -> Dict[Tuple[int, int], Form]:
        """
        Return EC curvature 2-forms :math:`R_{ab}` for *mode*.

        Returns
        -------
        dict mapping (a, b) → Form
        """
        C = self.structure_constants()
        Gamma_lc = self._koszul_connection(C)
        K = self._contortion(self.torsion_tensor_low(mode))
        Gamma = MutableDenseNDimArray.zeros(self.DIM, self.DIM, self.DIM)
        for a in range(self.DIM):
            for b in range(self.DIM):
                for c in range(self.DIM):
                    Gamma[a, b, c] = cancel(Gamma_lc[a, b, c] + K[a, b, c])
        R_abcd = self._lower_first_index(
            self._riemann_from_connection(Gamma, C)
        )
        forms: Dict[Tuple[int, int], Form] = {}
        for a in range(self.DIM):
            for b in range(self.DIM):
                form: Form = {}
                for c in range(self.DIM):
                    for d in range(self.DIM):
                        val = R_abcd[a, b, c, d] / 2
                        if val != S.Zero:
                            form = add_forms(
                                form,
                                scale_form(val, wedge(basis(c), basis(d))),
                            )
                forms[(a, b)] = form
        return forms

    def e_e_r_form(self, mode: str) -> Form:
        """Return :math:`e^a \\wedge e^b \\wedge R_{ab}` for *mode*."""
        R = self.curvature_forms_low(mode)
        out: Form = {}
        for a in range(self.DIM):
            for b in range(self.DIM):
                out = add_forms(
                    out, wedge(wedge(basis(a), basis(b)), R[(a, b)])
                )
        return out

    # ------------------------------------------------------------------
    # Full NY audit
    # ------------------------------------------------------------------

    def ny_form_audit(self, mode: str) -> dict:
        """
        Verify the Nieh-Yan exact-form identity for *mode*.

        Checks whether :math:`dQ_{NY} = T^a \\wedge T_a - e^a \\wedge e^b
        \\wedge R_{ab}` holds.

        Returns
        -------
        dict with keys:
          ``mode``, ``Q``, ``dQ``, ``TT``, ``e_e_R``,
          ``TT_minus_e_e_R``, ``TT_plus_e_e_R``,
          ``dQ_coeff``, ``TR_minus_coeff``, ``status``,
          ``integrated_Q``, ``reduced_dQ``.
        """
        qform = self.q_nieh_yan(mode)
        dqform = self.exterior_derivative(qform)
        tt = self.tt_form(mode)
        eer = self.e_e_r_form(mode)
        tr_minus = add_forms(tt, scale_form(-S.One, eer))
        tr_plus = add_forms(tt, eer)
        dq_coeff = coefficient(dqform)
        tr_minus_coeff = coefficient(tr_minus)

        if simplify(dq_coeff - tr_minus_coeff) == S.Zero:
            status = "NY_EXACT_CONFIRMED"
        elif simplify(dq_coeff + tr_minus_coeff) == S.Zero:
            status = "NY_EXACT_SIGN_FLIPPED"
        else:
            status = "NY_EXACT_CONVENTION_DEPENDENT"

        return {
            "mode": mode.upper(),
            "Q": qform,
            "dQ": dqform,
            "TT": tt,
            "e_e_R": eer,
            "TT_minus_e_e_R": tr_minus,
            "TT_plus_e_e_R": tr_plus,
            "dQ_coeff": dq_coeff,
            "TR_minus_coeff": tr_minus_coeff,
            "status": status,
            "integrated_Q": self.integrated_q_boundary(mode),
            "reduced_dQ": self.dq_reduced_density(mode),
        }
