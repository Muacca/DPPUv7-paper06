"""
Nil3Geometry — Nil³×S¹ Topology Implementation
================================================

EC+NY engine on Nil³ (Heisenberg nilmanifold) × S¹.

Characteristics of Nil³:
    - Only C²₀₁ is a nonzero structure constant  (1 axis vs. 6 axes for S³)
    - The radial coordinate follows the convention of using uppercase 'R'
      (an alias 'r' is also set for pipeline compatibility)
    - ε (squash) enters C²₀₁ as (1+ε)^{−4/3}
    - s (shear): volume-preserving shear of (e¹, e²) axes (enable_shear=True)
      Ansatz: e⁰=Rσ⁰, e¹=R(1+ε)^{2/3}(1+s)σ¹, e²=R(1+ε)^{−2/3}(1+s)^{−1}σ²
      → C²₀₁ = −(1+ε)^{−4/3}(1+s)^{−2}/R  (volume R³ preserved)
    - Under TWIST, since only dσ² = −σ⁰∧σ¹ is nonzero,
      **only ω₃** contributes to the structure constant C³₀₁ (ω₁, ω₂ are inactive)

Structure constants:
    Spatial block (diagonal): C²₀₁ = −(1+ε)^{−4/3}(1+s)^{−2} / R
    Spatial block (off-diag): C^a_{bc} via F = F_diag(ε,s) × G(z₃,z₄,z₅)  [enable_offdiag_shear]
    TWIST:                    C³₀₁ = −Lω₃ / R²

Volume factor: (2π)⁴ · L · R³   (volume-preserving: no ε or s factor)
"""

from sympy import pi, S, Rational
from sympy.tensor.array import MutableDenseNDimArray
from typing import Any, Dict

from .base_topology import TopologyEngine


class Nil3Geometry(TopologyEngine):
    """
    Engine dedicated to the Nil³×S¹ topology.

    Simpler than S³ but uses a different symbol naming convention ('R' vs 'r').

    Examples
    --------
    Direct usage::

        from dppu.topology import FiberMode, TopologyType, make_config
        from dppu.topology.nil3 import Nil3Geometry
        from dppu.torsion.mode import Mode

        cfg = make_config(
            TopologyType.NIL3,
            fiber_mode=FiberMode.TWIST,
            isotropic_twist=True,
            torsion_mode=Mode.AX,
        )
        engine = Nil3Geometry(cfg)
        engine.run()
    """

    # ── E4.1 hook ─────────────────────────────────────────────────────────────

    def _build_radial_and_deformation_params(self, params: Dict) -> None:
        """
        Add Nil³-specific symbols to params.

            R       : radial coordinate (positive — Nil³ convention)
            r       : alias for R (pipeline compatibility)
            epsilon : squash parameter (Symbol only when enable_squash=True)
            s       : shear parameter (Symbol when enable_shear=True; Zero otherwise)
                      Volume-preserving shear of (e¹, e²) axes:
                      e¹ = R(1+ε)^{2/3}(1+s) σ¹,  e² = R(1+ε)^{-2/3}(1+s)^{-1} σ²
                      → changes C²₀₁ but preserves volume R³
            q3,q4,q5: physical parameters for off-diagonal shear (for reference)
            z3,z4,z5: rationalized symbols z_i = exp(q_i/√2) when enable_offdiag_shear=True
                      Same Method-2 speedup as S³.  z=1 → G=I (diagonal frame).
        """
        from sympy import symbols, Symbol
        cfg = self.config

        R = symbols('R', positive=True, real=True)
        params['R'] = R
        params['r'] = R   # pipeline-compatibility alias for E4.4 (torsion ansatz), etc.

        params['epsilon'] = (
            symbols('epsilon', real=True) if cfg.enable_squash else S.Zero
        )
        params['s'] = (
            symbols('s', real=True) if cfg.enable_shear else S.Zero
        )

        enable_offdiag = getattr(cfg, 'enable_offdiag_shear', False)
        params['q3'] = symbols('q3', real=True) if enable_offdiag else S.Zero
        params['q4'] = symbols('q4', real=True) if enable_offdiag else S.Zero
        params['q5'] = symbols('q5', real=True) if enable_offdiag else S.Zero
        if enable_offdiag:
            params['z3'] = Symbol('z3', positive=True)
            params['z4'] = Symbol('z4', positive=True)
            params['z5'] = Symbol('z5', positive=True)
        else:
            params['z3'] = params['z4'] = params['z5'] = None

    # ── E4.2 hook ─────────────────────────────────────────────────────────────

    def _build_structure_constants(
        self, params: Dict, C: MutableDenseNDimArray
    ) -> None:
        """
        Write the Nil³ structure constants into C in-place.

        Fast path (q₃=q₄=q₅=0, diagonal frame):
            C²₀₁ = −(1+ε)^{−4/3}(1+s)^{−2} / R
            C²₁₀ = +(1+ε)^{−4/3}(1+s)^{−2} / R  (antisymmetric)

        General path (enable_offdiag_shear; F = F_diag × G):
            F_diag = diag(1, f₁, f₂), f₁=(1+ε)^{2/3}(1+s), f₂=(1+ε)^{-2/3}/(1+s)
            G = G₃(z₃)G₄(z₄)G₅(z₅)  [same hyperbolic-rotation ansatz as S³]
            Nil³ raw structure: only C^2_{01}(σ) = -1 nonzero
            → C^a_{bc}(e) = (−1/R)×F[a,2]×(F_inv[b,0]F_inv[c,1] − F_inv[b,1]F_inv[c,0])

        Additional TWIST terms (only ω₃ is active):
            C³₀₁ = −L ω₃ / R²
            C³₁₀ = +L ω₃ / R²
        """
        from .unified import FiberMode

        R       = params['R']
        L       = params['L']
        epsilon = params['epsilon']
        s       = params['s']

        # Spatial block
        z3 = params.get('z3')
        if z3 is None:
            # Fast path: diagonal frame only
            nil3_factor = (1 + epsilon) ** Rational(-4, 3) * (1 + s) ** (-2)
            C[2, 0, 1] = -nil3_factor / R
            C[2, 1, 0] = +nil3_factor / R
        else:
            # General path: F = F_diag_unit(ε,s) × G(z₃,z₄,z₅)
            # Method 2: z_i = exp(q_i/√2) — same rationalization as S³
            from sympy import Matrix, cancel as _cancel
            z4 = params['z4']
            z5 = params['z5']
            c3, s3_ = (z3 + 1/z3) / 2, (z3 - 1/z3) / 2
            c4, s4_ = (z4 + 1/z4) / 2, (z4 - 1/z4) / 2
            c5, s5_ = (z5 + 1/z5) / 2, (z5 - 1/z5) / 2

            G3 = Matrix([[c3,  s3_, S.Zero],
                         [s3_, c3,  S.Zero],
                         [S.Zero, S.Zero, S.One]])
            G4 = Matrix([[c4,  S.Zero, s4_],
                         [S.Zero, S.One, S.Zero],
                         [s4_, S.Zero, c4]])
            G5 = Matrix([[S.One, S.Zero, S.Zero],
                         [S.Zero, c5,  s5_],
                         [S.Zero, s5_, c5]])

            G = G3 * G4 * G5
            f1 = (1 + epsilon) ** Rational(2, 3)  * (1 + s)
            f2 = (1 + epsilon) ** Rational(-2, 3) / (1 + s)
            F_diag_unit = Matrix.diag(S.One, f1, f2)
            F = F_diag_unit * G

            G3_inv = Matrix([[c3,  -s3_, S.Zero],
                             [-s3_, c3,  S.Zero],
                             [S.Zero, S.Zero, S.One]])
            G4_inv = Matrix([[c4,  S.Zero, -s4_],
                             [S.Zero, S.One, S.Zero],
                             [-s4_, S.Zero, c4]])
            G5_inv = Matrix([[S.One, S.Zero, S.Zero],
                             [S.Zero, c5,  -s5_],
                             [S.Zero, -s5_, c5]])
            F_diag_unit_inv = Matrix.diag(S.One, S.One / f1, S.One / f2)
            F_inv = G5_inv * G4_inv * G3_inv * F_diag_unit_inv

            # Nil³ has only C^2_{01}(σ) = -1  →  2-term formula (vs. 6 terms for S³)
            # C^a_{bc}(e) = (-1/R) × F[a,2] × (F_inv[b,0]F_inv[c,1] - F_inv[b,1]F_inv[c,0])
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        v = F[a, 2] * (F_inv[b, 0] * F_inv[c, 1]
                                       - F_inv[b, 1] * F_inv[c, 0])
                        if v != S.Zero:
                            C[a, b, c] = _cancel(-v / R)

        # TWIST: In Nil³, only dσ² = −σ⁰∧σ¹ is nonzero → only ω₃ contributes
        if self.config.fiber_mode == FiberMode.TWIST:
            omega3 = params['omega3']  # iso=True → omega3=omega; iso=False → independent symbol
            twist_factor = L * omega3 / R ** 2
            C[3, 0, 1] = -twist_factor
            C[3, 1, 0] = +twist_factor

    def _compute_volume(self, params: Dict) -> Any:
        """(2π)⁴ · L · R³"""
        return (2 * pi) ** 4 * params['L'] * params['R'] ** 3
