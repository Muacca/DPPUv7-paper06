"""
T3Geometry — T³×S¹ Topology Implementation
============================================

EC+NY engine on T³ (3-dimensional torus, flat manifold) × S¹.

Characteristics:
    - All structure constants are zero  (C^a_{bc} = 0  for all a,b,c)
    - ε (squash): axisymmetric scale deformation  e^2 → r(1+ε) dx^2
                  Valid to introduce as a Symbol; structure constants remain exactly zero
                  (T³ is a flat manifold: Weyl curvature = 0 for any ε).
                  Used as a null test: SymPy confirms C²_LC = C²_EC = 0 ∀ε.
    - s (shear) : anisotropic deformation  e^2 → r(1+ε)(1+s) dx^2
                  Enabled by enable_shear=True; structure constants remain exactly zero.
                  Null test: SymPy confirms m²_shear = 0 for flat T³.
    - For TWIST, fiber symbols (ω₁,ω₂,ω₃) enter params but have no physical effect
      because dσ^i = 0, so they do not contribute to the structure constants
    - MIXING is specific to S³ and should be restricted at the DOFConfig level

Volume factor: (2π)⁴ · L · r³ · (1+ε) · (1+s)   [squash=True, shear=True]
               (2π)⁴ · L · r³ · (1+ε)             [squash=True, shear=False]
               (2π)⁴ · L · r³                      [squash=False, ε=s=0]
"""

from sympy import pi, S
from sympy.tensor.array import MutableDenseNDimArray
from typing import Any, Dict

from .base_topology import TopologyEngine


class T3Geometry(TopologyEngine):
    """
    Engine dedicated to the T³×S¹ topology.

    The simplest implementation: all structure constants are zero regardless of ε.
    When ``enable_squash=True``, an ``epsilon`` Symbol is introduced for the
    axisymmetric frame deformation  e^2 = r(1+ε) dx^2.  Because the structure
    constants are identically zero (T³ is flat), this serves as a null test:
    SymPy confirms C²_LC = C²_EC = 0 for arbitrary ε.

    Examples
    --------
    Direct usage::

        from dppu.topology.t3 import T3Geometry
        from dppu.topology import FiberMode, TopologyType, make_config
        from dppu.torsion.mode import Mode

        cfg = make_config(
            TopologyType.T3,
            fiber_mode=FiberMode.TWIST,
            isotropic_twist=False,
            torsion_mode=Mode.MX,
        )
        engine = T3Geometry(cfg)
        engine.run()
    """

    # ── E4.1 hook ─────────────────────────────────────────────────────────────

    def _build_radial_and_deformation_params(self, params: Dict) -> None:
        """
        Add T³-specific symbols to params.

            r       : radial coordinate (positive)
            epsilon : squash parameter (Symbol when enable_squash=True; Zero otherwise)
                      Represents the axisymmetric frame deformation  e^2 = r(1+ε) dx^2.
                      Does NOT appear in the structure constants (all remain zero).
                      Enters only the volume element as a (1+ε) factor.
            s       : shear parameter (Symbol when enable_shear=True; Zero otherwise)
                      Does NOT appear in structure constants (T³ flat: all C^a_{bc}=0).
                      Enters volume element as (1+s) factor.
                      Null test: ∂²V/∂s² = 0 confirms T³ spin-2 shear mass = 0.
        """
        from sympy import symbols
        cfg = self.config
        params['r'] = symbols('r', positive=True, real=True)
        params['epsilon'] = (
            symbols('epsilon', real=True) if cfg.enable_squash else S.Zero
        )
        params['s'] = (
            symbols('s', real=True) if cfg.enable_shear else S.Zero
        )

    # ── E4.2 hook ─────────────────────────────────────────────────────────────

    def _build_structure_constants(
        self, params: Dict, C: MutableDenseNDimArray
    ) -> None:
        """
        T³ is a flat manifold: all structure constants are zero for any ε or s.

        Neither squash (epsilon) nor shear (s) enter here; they only affect the
        volume element.  This verifies that C²_LC = C²_EC = 0 holds regardless
        of the metric deformation (null test for T³ flatness).
        """
        # T³: dσ^i = 0  →  C^a_{bc} = 0  for all indices, for any ε or s
        pass

    def _compute_volume(self, params: Dict) -> Any:
        """(2π)⁴ · L · r³ · (1+ε) · (1+s)

        The squash deformation  e^2 → r(1+ε) dx^2  and shear deformation
        (1+s) each scale the volume element independently.
        At ε=s=0 (or when enable_squash=enable_shear=False) this reduces to
        the isotropic result (2π)⁴ · L · r³.
        """
        return (2 * pi) ** 4 * params['L'] * params['r'] ** 3 * (1 + params['epsilon']) * (1 + params['s'])
