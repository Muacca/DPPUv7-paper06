"""
Base Frame Engine Pipeline
==========================

Abstract base class for Einstein-Cartan + Nieh-Yan computations.
Topology runners (S3S1, T3S1, Nil3S1) inherit from this.
"""

import traceback
from typing import Any, Dict, List, Optional, Tuple

from sympy import S, cancel, Rational, lambdify, Matrix
from sympy.tensor.array import MutableDenseNDimArray
import numpy as np

from ..torsion.mode import Mode
from ..torsion.nieh_yan import NyVariant
from ..utils.epsilon import epsilon_3d, epsilon_4d
from ..curvature.riemann import (
    verify_antisymmetry_strict, compute_riemann_tensor, lower_first_index
)
from ..utils.logger import NullLogger
from .checkpoint import CheckpointManager


class BaseFrameEngine:
    """
    Abstract Base Class for EC + NY Calculations.

    Topology-specific runners must implement:
    - step_E4_1_setup: Setup parameters
    - step_E4_2_metric_and_frame: Define metric, volume, structure constants
    """

    def __init__(
        self,
        mode: Mode,
        ny_variant: NyVariant,
        logger: Optional[Any] = None,
        checkpoint_mgr: Optional[CheckpointManager] = None,
        skip_antisymmetry_check: bool = False,
        weyl_source: Any = "ec",
        riemann_source: Any = "ec",
    ):
        if not isinstance(mode, Mode):
            raise TypeError(f"mode must be Mode enum, got {type(mode)}")
        if not isinstance(ny_variant, NyVariant):
            raise TypeError(f"ny_variant must be NyVariant enum")

        self.mode = mode
        self.ny_variant = ny_variant
        self.logger = logger or NullLogger()
        self.ckpt = checkpoint_mgr
        self.skip_antisymmetry_check = skip_antisymmetry_check
        self.weyl_source = self._normalize_curvature_source(weyl_source)
        self.riemann_source = self._normalize_curvature_source(riemann_source)
        self.data: Dict[str, Any] = {}

        self.steps = [
            ("E4.1", self.step_E4_1_setup),
            ("E4.2", self.step_E4_2_metric_and_frame),
            ("E4.3", self.step_E4_3_christoffel_frame),
            ("E4.3a", self.step_E4_3a_riemann_LC),
            ("E4.3b", self.step_E4_3b_weyl_tensor),
            ("E4.4", self.step_E4_4_torsion_ansatz_frame),
            ("E4.5", self.step_E4_5_contortion_frame),
            ("E4.6", self.step_E4_6_ec_connection_frame),
            ("E4.7", self.step_E4_7_riemann_tensor_frame),
            ("E4.8", self.step_E4_8_ricci_scalar_frame),
            ("E4.8b", self.step_E4_7b_weyl_ec),
            ("E4.9", self.step_E4_9_torsion_scalar_frame),
            ("E4.10", self.step_E4_10_nieh_yan_frame),
            ("E4.11", self.step_E4_11_lagrangian),
            ("E4.12", self.step_E4_12_angular_integration),
            ("E4.13", self.step_E4_13_effective_potential),
            ("E4.15", self.step_E4_15_summary),
        ]

    @staticmethod
    def _normalize_curvature_source(source: Any) -> str:
        """Normalize LC/EC curvature source enums or strings to lowercase keys."""
        value = getattr(source, "value", source)
        key = str(value).strip().lower()
        if key not in {"lc", "ec"}:
            raise ValueError(
                f"curvature source must be 'lc' or 'ec', got {source!r}"
            )
        return key

    def run(self, start_step: str = "E4.1"):
        """Execute computation pipeline."""
        self.logger.info(f"Mode: {self.mode.value}, NY Variant: {self.ny_variant.value}")

        try:
            start_idx = next(
                i for i, (sid, _) in enumerate(self.steps) if sid == start_step
            )
        except StopIteration:
            raise ValueError(f"Invalid start step: {start_step}")

        for step_id, step_func in self.steps[start_idx:]:
            try:
                doc = step_func.__doc__ or "Executing step..."
                self.logger.step(step_id, doc.strip())
                step_func()
                if self.ckpt:
                    self.ckpt.save(step_id, self.data)
            except Exception as e:
                self.logger.error(f"Step {step_id} failed: {e}")
                raise

        self.logger.finalize()

    # Abstract methods
    def step_E4_1_setup(self):
        """Setup parameters - MUST be implemented by topology runner."""
        raise NotImplementedError

    def step_E4_2_metric_and_frame(self):
        """Define metric/volume/structure constants - MUST be implemented."""
        raise NotImplementedError

    def step_E4_3_christoffel_frame(self):
        """Compute LC connection via Koszul formula."""
        from .levi_civita import compute_christoffel_frame
        from .metric import verify_metric_compatibility

        dim = self.data['dim']
        C = self.data['structure_constants']

        Gamma_LC = compute_christoffel_frame(C, dim, self.logger)
        self.data['connection_LC'] = Gamma_LC

        # Verify metric compatibility
        verify_metric_compatibility(
            Gamma_LC, self.data['metric_frame'], dim, self.logger
        )

    def step_E4_3a_riemann_LC(self):
        """Compute Riemann/Ricci tensors for Levi-Civita connection."""
        from ..curvature.ricci import compute_ricci_tensor, compute_ricci_scalar

        dim = self.data['dim']
        Gamma_LC = self.data['connection_LC']
        C = self.data['structure_constants']

        self.logger.info("Computing LC Riemann Tensor...")
        # Note: compute_riemann_tensor uses generic Gamma. Pass Gamma_LC here.
        Riemann_LC = compute_riemann_tensor(Gamma_LC, C, dim)
        self.data['riemann_LC'] = Riemann_LC

        # Lower indices and verify antisymmetry
        R_abcd_LC = lower_first_index(Riemann_LC, self.data['metric_frame'], dim)
        if not self.skip_antisymmetry_check:
            verify_antisymmetry_strict(R_abcd_LC, dim, None, self.logger)
        else:
            self.logger.info("LC Riemann antisymmetry check: SKIPPED")
        self.data['riemann_abcd_LC'] = R_abcd_LC

        # Compute Ricci components for LC (needed for Weyl)
        self.data['ricci_LC'] = compute_ricci_tensor(
            Riemann_LC, dim, self.logger
        )
        self.data['ricci_scalar_LC'] = compute_ricci_scalar(
            Riemann_LC, dim, self.logger
        )

    def step_E4_3b_weyl_tensor(self):
        """Compute Weyl tensor and scalar."""
        from ..curvature.weyl import compute_weyl_tensor, compute_weyl_scalar

        dim = self.data['dim']
        
        # Inputs from LC step
        R_abcd = self.data['riemann_abcd_LC']
        Ricci = self.data['ricci_LC']
        R_scalar = self.data['ricci_scalar_LC']
        metric = self.data['metric_frame']
        
        # Compute Weyl Tensor
        C_abcd = compute_weyl_tensor(
            R_abcd, Ricci, R_scalar, metric, dim, self.logger
        )
        self.data['weyl_tensor'] = C_abcd
        
        # Compute Weyl Scalar
        # Note: frame metric inverse is same as metric if diag(1,1,1,1)
        # But generally we should invoke inverse.
        metric_inv = metric.inv() 
        C_sq = compute_weyl_scalar(C_abcd, metric_inv, dim, self.logger)
        self.data['weyl_scalar'] = C_sq

    def step_E4_4_torsion_ansatz_frame(self):
        """Construct torsion tensor."""
        from ..torsion.ansatz import construct_torsion_tensor

        dim = self.data['dim']
        r = self.data['params']['r']
        eta = self.data['params']['eta']
        V = self.data['params']['V']
        metric = self.data['metric_frame']
        fc  = getattr(self.config, 'frame_convention', 'legacy_euclidean')
        sig = getattr(self.config, 'signature', 'euclidean')

        T = construct_torsion_tensor(
            self.mode, r, eta, V, metric, dim, self.logger,
            frame_convention=fc,
            signature=sig,
        )
        self.data['torsion_tensor'] = T

    def step_E4_5_contortion_frame(self):
        """Compute contortion from torsion."""
        from .contortion import compute_contortion

        dim    = self.data['dim']
        T      = self.data['torsion_tensor']
        metric = self.data['metric_frame']
        K = compute_contortion(T, dim, metric=metric, logger=self.logger)
        self.data['contortion'] = K

    def step_E4_6_ec_connection_frame(self):
        """Compute EC connection = LC + contortion."""
        from .ec_connection import compute_ec_connection

        dim = self.data['dim']
        Gamma_LC = self.data['connection_LC']
        K = self.data['contortion']
        Gamma_EC = compute_ec_connection(Gamma_LC, K, dim, self.logger)
        self.data['connection_EC'] = Gamma_EC

    def step_E4_7_riemann_tensor_frame(self):
        """Compute Riemann tensor with antisymmetry verification."""
        dim = self.data['dim']
        Gamma = self.data['connection_EC']
        C = self.data['structure_constants']

        self.logger.info("Computing EC Riemann Tensor...")
        Riemann = compute_riemann_tensor(Gamma, C, dim)
        self.data['riemann'] = Riemann

        # Lower first index and verify antisymmetry
        R_abcd = lower_first_index(Riemann, self.data['metric_frame'], dim)
        if not self.skip_antisymmetry_check:
            verify_antisymmetry_strict(R_abcd, dim, None, self.logger)
        else:
            self.logger.info("EC Riemann antisymmetry check: SKIPPED")
        self.data['riemann_abcd'] = R_abcd

        self.logger.success("Riemann Tensor computed and verified")

    def step_E4_7b_weyl_ec(self):
        """Compute Weyl scalar for EC connection (C²_EC) and store as weyl_scalar_EC.

        Requires: data['riemann'], data['riemann_abcd'], data['ricci_scalar'],
                  data['metric_frame'], data['dim']  (all available after E4.8).

        C²_EC = C²_LC + 16V²η²/(3r²) for MX mode; equals C²_LC for AX/VT (Theorem 1).
        data['weyl_scalar'] (C²_LC, from E4.3b) is kept unchanged so that proof
        scripts can still compare LC vs EC Weyl scalars explicitly.
        """
        from ..action.ec_action import compute_c2_ec
        self.data['weyl_scalar_EC'] = compute_c2_ec(self.data)

    def step_E4_8_ricci_scalar_frame(self):
        """Compute Ricci scalar."""
        from ..curvature.ricci import compute_ricci_scalar

        dim = self.data['dim']
        R_scalar = compute_ricci_scalar(self.data['riemann'], dim, self.logger)
        self.data['ricci_scalar'] = R_scalar

    def step_E4_9_torsion_scalar_frame(self):
        """Compute torsion scalar."""
        from ..torsion.scalar import compute_torsion_scalar

        dim = self.data['dim']
        T_scalar = compute_torsion_scalar(
            self.data['torsion_tensor'], dim, self.logger
        )
        self.data['torsion_scalar'] = T_scalar

    def step_E4_10_nieh_yan_frame(self):
        """Compute Nieh-Yan density (all variants)."""
        dim = self.data['dim']
        T = self.data['torsion_tensor']
        R_abcd = self.data['riemann_abcd']

        self.logger.info(f"Computing Nieh-Yan (variant: {self.ny_variant.value})...")

        # TT term
        N_TT = S.Zero
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):
                        eps = epsilon_4d(a, b, c, d)
                        if eps != 0:
                            for e in range(dim):
                                N_TT += Rational(1, 4) * eps * T[e, a, b] * T[e, c, d]
        self.data['ny_density_TT'] = cancel(N_TT)

        # Ree term
        N_Ree = S.Zero
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):
                        eps = epsilon_4d(a, b, c, d)
                        if eps != 0:
                            N_Ree += Rational(1, 4) * eps * R_abcd[a, b, c, d]
        self.data['ny_density_Ree'] = cancel(N_Ree)

        # Full
        self.data['ny_density_full'] = cancel(N_TT - N_Ree)

        # Select adopted variant
        if self.ny_variant == NyVariant.TT:
            self.data['nieh_yan_density'] = self.data['ny_density_TT']
        elif self.ny_variant == NyVariant.REE:
            self.data['nieh_yan_density'] = self.data['ny_density_Ree']
        else:  # FULL
            self.data['nieh_yan_density'] = self.data['ny_density_full']

        self.logger.success(f"Nieh-Yan: {self.ny_variant.value}")

    def step_E4_11_lagrangian(self):
        """Construct Lagrangian."""
        from ..action.lagrangian import compute_lagrangian

        weyl_key = 'weyl_scalar' if self.weyl_source == 'lc' else 'weyl_scalar_EC'
        weyl_scalar = self.data.get(weyl_key, S.Zero)
        self.data['weyl_scalar_action'] = weyl_scalar
        self.data['weyl_source_action'] = self.weyl_source
        self.logger.info(f"Weyl source for action: {self.weyl_source.upper()}")

        L = compute_lagrangian(
            self.data['ricci_scalar'],
            self.data['nieh_yan_density'],
            self.data['params']['kappa'],
            self.data['params']['theta_NY'],
            weyl_scalar=weyl_scalar,
            alpha=self.data['params'].get('alpha', S.Zero),
            logger=self.logger
        )
        self.data['lagrangian'] = L

    def step_E4_12_angular_integration(self):
        """Integrate action."""
        from ..action.lagrangian import compute_action

        action = compute_action(
            self.data['lagrangian'],
            self.data['total_volume'],
            self.logger
        )
        self.data['action'] = action

    def step_E4_13_effective_potential(self):
        """Extract effective potential."""
        from ..action.potential import compute_effective_potential

        V = compute_effective_potential(self.data['action'], self.logger)
        self.data['potential'] = V

    def step_E4_15_summary(self):
        """Summary."""
        self.logger.info("Computation complete.")

    def get_effective_potential_function(self):
        """Get lambdified V(r) for numerical evaluation.

        Returns a callable with signature:
            f(r, V, eta, theta_NY, L, kappa, epsilon, alpha) -> float

        Parameters that are fixed constants (e.g. epsilon=0 for T3) are
        accepted but ignored, so callers always use the same 8-arg signature.
        """
        if 'potential' not in self.data:
            raise RuntimeError("Run engine first.")

        from sympy import Symbol
        params = self.data['params']
        _KEYS = ['r', 'V', 'eta', 'theta_NY', 'L', 'kappa', 'epsilon', 'alpha']
        all_vals = [params[k] for k in _KEYS]

        # Only Symbol entries become lambdify arguments
        active_indices = [i for i, v in enumerate(all_vals) if isinstance(v, Symbol)]
        active_args = [all_vals[i] for i in active_indices]

        inner_func = lambdify(active_args, self.data['potential'], modules='numpy', cse=True)

        def wrapper(r, V, eta, theta_NY, L, kappa, epsilon, alpha):
            full = [r, V, eta, theta_NY, L, kappa, epsilon, alpha]
            return inner_func(*[full[i] for i in active_indices])

        return wrapper

    def get_potential_decomposition(self) -> Dict:
        """Decompose potential by powers of r."""
        from ..action.potential import decompose_potential
        return decompose_potential(self.data['potential'], self.data['params']['r'])
