"""
Sol3Geometry -- Sol x S1 Topology Implementation
=================================================

EC+NY engine on Sol (solvable Lie group, one of Thurston's 8 geometries) x S1.

Characteristics of Sol:
    - Two nonzero structure constants in the spatial block:
          C[1, 0, 1] = +1/R   (from de^1 = e^0 ^ e^1 / R)
          C[2, 0, 2] = -1/R   (from de^2 = -e^0 ^ e^2 / R)
      (and antisymmetric partners)
    - KEY PROPERTY: These are INDEPENDENT of frame scaling (epsilon, s).
      The exponential factors e^z, e^{-z} in the Sol coframe cancel the
      scaling, making structure constants rigid functions of R only.
    - Consequence: Riemann tensor, Weyl scalar are epsilon, s-independent.
      For FiberMode.NONE: epsilon, s do NOT enter V_eff at all
      (volume is also rigid: f1*f2 = 1, confirmed by SymPy auto-simplification).
      For FiberMode.TWIST: epsilon, s enter through C[3,0,1] and C[3,0,2] only.

Lie algebra (DPPU convention: [Eb, Ec] = -C^a_{bc} Ea):
    [E0, E1] = -E1/R  ->  C^1_{01} = +1/R
    [E0, E2] = +E2/R  ->  C^2_{02} = -1/R
    (All other spatial commutators vanish)

Coframe (unit Sol, R=1):
    sigma^0 = dz
    sigma^1 = e^z dx
    sigma^2 = e^{-z} dy

    d sigma^0 = 0                  -> C[0, ...] = 0 in spatial block
    d sigma^1 = sigma^0 ^ sigma^1  -> C[1, 0, 1] = +1/R
    d sigma^2 = -sigma^0 ^ sigma^2 -> C[2, 0, 2] = -1/R

Deformed coframe (frame scaling epsilon, s, same form as Nil3 for fair comparison):
    e^0 = R sigma^0
    e^1 = R (1+eps)^{2/3} (1+s) sigma^1   [f1 = (1+eps)^{2/3}(1+s)]
    e^2 = R (1+eps)^{-2/3} / (1+s) sigma^2 [f2 = (1+eps)^{-2/3}/(1+s)]

    Because e^z in sigma^1 cancels f1:
        de^1 = R*f1 * d sigma^1 = R*f1 * sigma^0^sigma^1
             = R*f1 * (e^0/R) ^ (e^1/(R*f1)) = e^0^e^1/R  (f1-independent!)
    Similarly de^2 = -e^0^e^2/R  (f2-independent)
    -> Structure constants remain C[1,0,1]=+1/R, C[2,0,2]=-1/R for any eps, s.

Volume factor: (2*pi)^4 * L * R^3   (volume-preserving in eps, s; same as Nil3)

TWIST fiber terms (when FiberMode.TWIST):
    sigma^1 = e^z dx, sigma^2 = e^{-z} dy, so:
        d sigma^1 = sigma^0 ^ sigma^1   -> omega1 contributes to C[3, 0, 1]
                                           Wait: only dσ^1, dσ^2 are nonzero.
    For Sol TWIST: e^3 = L(dt + w1*sigma^0 + w2*sigma^1 + w3*sigma^2)
        de^3 = L*(w1*d sigma^0 + w2*d sigma^1 + w3*d sigma^2)
             = L*(0 + w2*sigma^0^sigma^1 + w3*(-sigma^0^sigma^2))
             = L*w2*(e^0/R)^(e^1/(R*f1)) - L*w3*(e^0/R)^(e^2/(R*f2))
             = L*w2/(R^2*f1) * e^0^e^1 - L*w3/(R^2*f2) * e^0^e^2
    -> C[3, 0, 1] = L*w2 / (R^2*f1)   [depends on eps, s through f1]
       C[3, 0, 2] = -L*w3 / (R^2*f2)  [depends on eps, s through f2]
    (Unlike spatial block, TWIST terms ARE eps, s-dependent for Sol.)
    Note: w1 does not contribute since d sigma^0 = 0.
"""

from sympy import cancel, pi, S, Rational
from sympy.tensor.array import MutableDenseNDimArray
from typing import Any, Dict

from .base_topology import TopologyEngine


class Sol3Geometry(TopologyEngine):
    """
    Engine dedicated to the Sol x S1 topology.

    Sol is the unique solvable Lie geometry among Thurston's 8 geometries
    that has no compact quotient of the standard 3-sphere type. Its key
    structural property in this pipeline is that the spatial structure
    constants C[1,0,1]=+1/R and C[2,0,2]=-1/R are INDEPENDENT of frame
    deformation parameters (epsilon, s), unlike S3 or Nil3.

    Examples
    --------
    Direct usage::

        from dppu.topology.sol3 import Sol3Geometry
        from dppu.topology.unified import DOFConfig, TopologyType
        from dppu.torsion.mode import Mode

        cfg = DOFConfig(
            topology=TopologyType.SOL3,
            torsion_mode=Mode.AX,
            enable_squash=True,
            enable_shear=True,
        )
        engine = Sol3Geometry(cfg)
        engine.run()
    """

    # -- E4.1 hook -----------------------------------------------------------

    def _build_radial_and_deformation_params(self, params: Dict) -> None:
        """
        Add Sol-specific symbols to params.

            R       : radial coordinate (positive -- Sol convention, same as Nil3)
            r       : alias for R (pipeline compatibility)
            epsilon : squash parameter (Symbol when enable_squash=True; Zero otherwise)
                      NOTE: does NOT affect structure constants (Sol rigidity property).
                      Enters only the torsion ansatz via r=R and potentially the
                      TWIST fiber terms.
            s       : shear parameter (Symbol when enable_shear=True; Zero otherwise)
                      NOTE: does NOT affect spatial structure constants (Sol rigidity).
        """
        from sympy import symbols
        cfg = self.config

        R = symbols('R', positive=True, real=True)
        params['R'] = R
        params['r'] = R   # pipeline-compatibility alias

        params['epsilon'] = (
            symbols('epsilon', real=True) if cfg.enable_squash else S.Zero
        )
        params['s'] = (
            symbols('s', real=True) if cfg.enable_shear else S.Zero
        )

    # -- E4.2 hook -----------------------------------------------------------

    def _build_structure_constants(
        self, params: Dict, C: MutableDenseNDimArray
    ) -> None:
        """
        First-principles computation of Sol structure constants.

        For a diagonal frame deformation e^a = lambda_a * sigma^a, the structure
        constants transform as:
            C[a, b, c] = (lambda_a / (lambda_b * lambda_c)) * C_unit[a, b, c]
        where C_unit are unit-Sol (R=1, eps=s=0) structure constants:
            C_unit[1, 0, 1] = +1   (from d sigma^1 = sigma^0 ^ sigma^1)
            C_unit[2, 0, 2] = -1   (from d sigma^2 = -sigma^0 ^ sigma^2)

        Frame factors for deformed Sol:
            lambda_0 = R,  lambda_1 = R*f1,  lambda_2 = R*f2

        Spatial block (SymPy cancel() computes the cancellation explicitly):
            C[1,0,1] = (R*f1) / (R * R*f1) = 1/R   [f1 cancels: Sol rigidity]
            C[2,0,2] = -(R*f2) / (R * R*f2) = -1/R  [f2 cancels: Sol rigidity]

        TWIST fiber terms (lambda_3 = L; f1,f2 do NOT cancel here):
            C[3, 0, 1] = L*omega2 / (R * R*f1)
            C[3, 0, 2] = -L*omega3 / (R * R*f2)
        """
        from .unified import FiberMode

        R       = params['R']
        epsilon = params['epsilon']
        s       = params['s']

        # Frame scaling factors for the deformed Sol coframe:
        #   e^0 = R * sigma^0              (lambda_0 = R)
        #   e^1 = R * f1 * sigma^1         (lambda_1 = R * f1)
        #   e^2 = R * f2 * sigma^2         (lambda_2 = R * f2)
        f1 = (1 + epsilon) ** Rational(2, 3) * (1 + s)
        f2 = (1 + epsilon) ** Rational(-2, 3) / (1 + s)
        lambda_0 = R
        lambda_1 = R * f1
        lambda_2 = R * f2

        # Spatial block: C[a,b,c] = (lambda_a / (lambda_b * lambda_c)) * C_unit[a,b,c]
        # C[1,0,1] = (lambda_1 / (lambda_0 * lambda_1)) * (+1) = 1/R  [f1 cancels]
        val_101 = cancel(lambda_1 / (lambda_0 * lambda_1))
        C[1, 0, 1] = +val_101
        C[1, 1, 0] = -val_101   # antisymmetry

        # C[2,0,2] = (lambda_2 / (lambda_0 * lambda_2)) * (-1) = -1/R  [f2 cancels]
        val_202 = cancel(-lambda_2 / (lambda_0 * lambda_2))
        C[2, 0, 2] = +val_202
        C[2, 2, 0] = -val_202   # antisymmetry

        # TWIST fiber terms: omega2 -> C[3,0,1], omega3 -> C[3,0,2]
        # (omega1 inactive because d sigma^0 = 0 for Sol)
        if self.config.fiber_mode == FiberMode.TWIST:
            L      = params['L']
            omega2 = params['omega2']
            omega3 = params['omega3']

            # lambda_3 = L; C[3,0,i] = (L / (lambda_0 * lambda_i)) * omega_i
            # f1, f2 in lambda_1, lambda_2 do NOT cancel (no Sol rigidity for TWIST)
            if omega2 != S.Zero:
                twist_01 = L * omega2 / (lambda_0 * lambda_1)
                C[3, 0, 1] = +twist_01
                C[3, 1, 0] = -twist_01

            if omega3 != S.Zero:
                twist_02 = L * omega3 / (lambda_0 * lambda_2)
                C[3, 0, 2] = -twist_02
                C[3, 2, 0] = +twist_02

    def _compute_volume(self, params: Dict) -> Any:
        """Physical volume from first principles: det(frame) * (2*pi)^4 * L.

        V = (2*pi)^4 * L * lambda_0 * lambda_1 * lambda_2
          = (2*pi)^4 * L * R * (R*f1) * (R*f2)
          = (2*pi)^4 * L * R^3 * f1*f2

        f1 * f2 = (1+eps)^{2/3-2/3} * (1+s)^{1-1} = 1 (SymPy auto-simplifies),
        so the result is (2*pi)^4 * L * R^3 -- volume-preserving confirmed.
        """
        L = params['L']
        R = params['R']
        epsilon = params['epsilon']
        s = params['s']

        f1 = (1 + epsilon) ** Rational(2, 3) * (1 + s)
        f2 = (1 + epsilon) ** Rational(-2, 3) / (1 + s)

        lambda_0 = R
        lambda_1 = R * f1
        lambda_2 = R * f2

        # f1 * f2 = 1 (SymPy auto-simplifies), so this returns (2*pi)^4 * L * R^3
        return (2 * pi) ** 4 * L * lambda_0 * lambda_1 * lambda_2
