"""
TopologyEngine — Abstract Base for Topology-Specific EC+NY Engines
===================================================================

Separates the "common processing" and "topology-specific processing"
of E4.1/E4.2 using the Template Method pattern.

The full pipeline (E4.3–E4.15) is provided as-is by BaseFrameEngine.

Inheritance hierarchy
---------------------
    BaseFrameEngine          (dppu/engine/pipeline.py)
        └── TopologyEngine   (this class)  ← abstract
                ├── S3Geometry    (s3.py)
                ├── T3Geometry    (t3.py)
                └── Nil3Geometry  (nil3.py)

3 abstract methods that subclasses must implement
-------------------------------------------------
    _build_radial_and_deformation_params(params)
        → Add r/R, epsilon, s to params.

    _build_structure_constants(params, C)
        → Fill the structure constant matrix C in-place,
          including spatial block + fiber terms.

    _compute_volume(params) -> sympy expression
        → Return the volume factor arising from the angular integral
          over the 4-dimensional manifold.

Common methods (subclasses normally do not need to override)
------------------------------------------------------------
    _build_common_params(params)       L, kappa, theta_NY, alpha
    _build_fiber_params(params)        omega1-3 / delta
    _build_torsion_params(params)      eta, V
    _add_s3_twist_C(C, ...)            Embed S³ twist terms in C
    _add_s3_mixing_C(C, ...)           Embed S³ mixing terms in C
    get_riemann_lambdified()           lambdify for numerical Pontryagin scanning
    get_free_params()                  dict of free symbols
    get_config_summary()               one-line configuration summary string
"""

from abc import abstractmethod
from typing import Any, Dict, Optional

from sympy import symbols, Matrix, S, pi, Rational, sqrt
from sympy.tensor.array import MutableDenseNDimArray

from ..engine.pipeline import BaseFrameEngine
from ..torsion.mode import Mode
from ..torsion.nieh_yan import NyVariant
from ..utils.epsilon import epsilon_3d


class TopologyEngine(BaseFrameEngine):
    """
    Abstract base for topology-specific EC+NY engines.

    Constructs E4.1/E4.2 using the Template Method pattern.
    The remaining 13 pipeline steps (E4.3–E4.15) are provided by BaseFrameEngine.
    """

    def __init__(
        self,
        config,                          # DOFConfig  (type annotation omitted to avoid circular import)
        logger: Optional[Any] = None,
        checkpoint_mgr=None,
    ):
        # DOFConfig type check is performed in UnifiedEngine.__new__, so omitted here
        self.config = config
        super().__init__(
            mode=config.torsion_mode,
            ny_variant=config.ny_variant,
            logger=logger,
            checkpoint_mgr=checkpoint_mgr,
            skip_antisymmetry_check=getattr(config, 'skip_antisymmetry_check', False),
            weyl_source=getattr(config, 'weyl_source', 'ec'),
            riemann_source=getattr(config, 'riemann_source', 'ec'),
        )

    # =========================================================================
    # E4.1 — Parameter setup  (template method)
    # =========================================================================

    def step_E4_1_setup(self):
        """
        Define SymPy symbols for active DOFs; disabled DOFs are set to Zero.

        Processing order:
          1. _build_common_params        common symbols (L, kappa, theta_NY, alpha)
          2. _build_radial_and_deformation_params   ★ abstract  (r/R, ε, s)
          3. _build_fiber_params         fiber symbols (ω₁ω₂ω₃ or δ or Zero)
          4. _build_torsion_params       torsion symbols (η, V)
        """
        params: Dict[str, Any] = {}

        self._build_common_params(params)
        self._build_radial_and_deformation_params(params)
        self._build_fiber_params(params)
        self._build_torsion_params(params)

        params['q'] = 2 * params['eta']   # Keep pipeline compatibility

        self.data['params'] = params
        self.data['dim']    = 4

        self._log_setup_summary()

    # =========================================================================
    # E4.2 — Metric, volume, structure constants  (template method)
    # =========================================================================

    def step_E4_2_metric_and_frame(self):
        """
        Build metric, volume element, and structure constants.

        Processing order:
          1. Build frame metric (Euclidean identity or Lorentzian diag(-,+,+,+))
          2. _build_structure_constants   ★ abstract  (C in-place)
          3. _compute_volume              ★ abstract  → total_volume
          4. Store metric metadata in engine.data
        """
        from ..engine.metric import create_frame_metric
        from sympy import sign as sympy_sign

        params = self.data['params']
        dim    = self.data['dim']

        sig = getattr(self.config, 'signature', 'euclidean')
        fc  = getattr(self.config, 'frame_convention', 'legacy_euclidean')

        metric = create_frame_metric(dim, signature=sig)
        self.data['metric_frame']          = metric
        self.data['metric_inv_frame']      = metric.inv()
        self.data['metric_det']            = metric.det()
        self.data['metric_signature_sign'] = sympy_sign(metric.det())
        self.data['signature']             = sig
        self.data['frame_convention']      = fc
        self.data['volume_convention'] = (
            "legacy_euclidean" if sig == "euclidean"
            else "lorentzian_smoke_static_not_physical_action"
        )

        C = MutableDenseNDimArray.zeros(dim, dim, dim)
        self._build_structure_constants(params, C)

        if sig == "lorentzian" and fc == "lz_native":
            from .lz_adapter import remap_dense_structure_constants_to_lz
            C = remap_dense_structure_constants_to_lz(C, dim)

        self.data['structure_constants'] = C

        self.data['total_volume'] = self._compute_volume(params)

        self.logger.success(
            f"{self.config.topology.value.upper()}×S¹: "
            f"structure constants defined  "
            f"[fiber={self.config.fiber_mode.value}]  "
            f"[sig={sig}]"
        )

    # =========================================================================
    # Abstract hooks  (must be implemented in subclasses)
    # =========================================================================

    @abstractmethod
    def _build_radial_and_deformation_params(self, params: Dict) -> None:
        """
        Add radial symbols and deformation parameters to params.

        Implementation example (S³):
            params['r']       = symbols('r', positive=True)
            params['epsilon'] = symbols('epsilon') if enable_squash else S.Zero
            params['s']       = symbols('s')       if enable_shear  else S.Zero
        """

    @abstractmethod
    def _build_structure_constants(
        self, params: Dict, C: MutableDenseNDimArray
    ) -> None:
        """
        Fill the structure constant matrix C in-place.

        Includes the spatial block (i,j,k ∈ {0,1,2}) and fiber terms
        (TWIST: C³_{jk}, MIXING: C⁰_{i3}, etc.).
        For T³, nothing needs to be done (C remains zero).
        """

    @abstractmethod
    def _compute_volume(self, params: Dict) -> Any:
        """
        Return the SymPy expression for the volume factor arising from angular integration.

        Examples:
            S³:   2 * pi**2 * L * r**3
            T³:   (2*pi)**4 * L * r**3
            Nil³: (2*pi)**4 * L * R**3
        """

    # =========================================================================
    # Concrete helpers — common to all topologies
    # =========================================================================

    def _build_common_params(self, params: Dict) -> None:
        """Add the 4 symbols that are always present for every topology."""
        params['L']        = symbols('L',        positive=True, real=True)
        params['kappa']    = symbols('kappa',    positive=True, real=True)
        params['theta_NY'] = symbols('theta_NY', real=True)
        params['alpha']    = symbols('alpha',    real=True)

    def _build_fiber_params(self, params: Dict) -> None:
        """
        Add fiber symbols according to FiberMode.

        TWIST (iso)    → omega (= omega1 = omega2 = omega3); delta0=delta1=delta2=0
        TWIST (aniso)  → omega1, omega2, omega3;             delta0=delta1=delta2=0
        MIXING         → delta0, delta1, delta2;             omega1=omega2=omega3=0
        BOTH           → delta0, delta1, delta2 + omega1, omega2, omega3
        NONE           → all Zero

        Note: delta is fully vectorized to a 3-component vector (δ₀,δ₁,δ₂).
              The old scalar 'delta' symbol is deprecated. Use delta2 in scripts.
        """
        # Avoid circular import: FiberMode is imported from unified
        from .unified import FiberMode
        cfg = self.config

        if cfg.fiber_mode == FiberMode.TWIST:
            if cfg.isotropic_twist:
                omega = symbols('omega', real=True)
                params['omega']  = omega
                params['omega1'] = omega
                params['omega2'] = omega
                params['omega3'] = omega
            else:
                params['omega1'] = symbols('omega1', real=True)
                params['omega2'] = symbols('omega2', real=True)
                params['omega3'] = symbols('omega3', real=True)
            params['delta0'] = S.Zero
            params['delta1'] = S.Zero
            params['delta2'] = S.Zero

        elif cfg.fiber_mode == FiberMode.MIXING:
            # 3-component vector δ = (δ₀,δ₁,δ₂)
            # δ₀: rotation in the (ê⁰,ê³) plane, δ₁: (ê¹,ê³) plane, δ₂: (ê²,ê³) plane
            params['delta0'] = symbols('delta0', real=True)
            params['delta1'] = symbols('delta1', real=True)
            params['delta2'] = symbols('delta2', real=True)
            params['omega1'] = S.Zero
            params['omega2'] = S.Zero
            params['omega3'] = S.Zero
            # Polynomial auxiliary symbols: cd_i = cos θ_i = 1/√(1+δ_i²), sd_i = sin θ_i = δ_i/√(1+δ_i²)
            # Keep polynomial arithmetic throughout the pipeline to avoid rational expressions
            params['cd0'] = symbols('cd0', real=True)
            params['sd0'] = symbols('sd0', real=True)
            params['cd1'] = symbols('cd1', real=True)
            params['sd1'] = symbols('sd1', real=True)
            params['cd2'] = symbols('cd2', real=True)
            params['sd2'] = symbols('sd2', real=True)

        elif cfg.fiber_mode == FiberMode.BOTH:
            # TWIST + MIXING simultaneously. isotropic_twist is unsupported (always anisotropic).
            params['delta0'] = symbols('delta0', real=True)
            params['delta1'] = symbols('delta1', real=True)
            params['delta2'] = symbols('delta2', real=True)
            params['omega1'] = symbols('omega1', real=True)
            params['omega2'] = symbols('omega2', real=True)
            params['omega3'] = symbols('omega3', real=True)
            # Polynomial auxiliary symbols: cd_i = cos θ_i = 1/√(1+δ_i²), sd_i = sin θ_i = δ_i/√(1+δ_i²)
            params['cd0'] = symbols('cd0', real=True)
            params['sd0'] = symbols('sd0', real=True)
            params['cd1'] = symbols('cd1', real=True)
            params['sd1'] = symbols('sd1', real=True)
            params['cd2'] = symbols('cd2', real=True)
            params['sd2'] = symbols('sd2', real=True)
            # Velocity symbols (for G metric computation; Symbol only when enable_velocity=True)
            if getattr(cfg, 'enable_velocity', False):
                params['v_omega1'] = symbols('v_omega1', real=True)
                params['v_omega2'] = symbols('v_omega2', real=True)
                params['v_omega3'] = symbols('v_omega3', real=True)
                params['v_delta0'] = symbols('v_delta0', real=True)
                params['v_delta1'] = symbols('v_delta1', real=True)
                params['v_delta2'] = symbols('v_delta2', real=True)
            else:
                params['v_omega1'] = S.Zero
                params['v_omega2'] = S.Zero
                params['v_omega3'] = S.Zero
                params['v_delta0'] = S.Zero
                params['v_delta1'] = S.Zero
                params['v_delta2'] = S.Zero

        else:  # FiberMode.NONE
            params['delta0'] = S.Zero
            params['delta1'] = S.Zero
            params['delta2'] = S.Zero
            params['omega1'] = S.Zero
            params['omega2'] = S.Zero
            params['omega3'] = S.Zero

    def _build_torsion_params(self, params: Dict) -> None:
        """Add eta, V symbols according to Mode."""
        if self.mode == Mode.AX:
            params['eta'] = symbols('eta', real=True)
            params['V']   = S.Zero
        elif self.mode == Mode.VT:
            params['eta'] = S.Zero
            params['V']   = symbols('V', positive=True, real=True)
        else:  # MX
            params['eta'] = symbols('eta', real=True)
            params['V']   = symbols('V',   positive=True, real=True)

    def _log_setup_summary(self) -> None:
        """Log the list of active DOFs."""
        from .unified import FiberMode
        cfg = self.config
        active = []
        if cfg.enable_squash:                     active.append("ε")
        if cfg.enable_shear:                      active.append("s")
        if cfg.fiber_mode == FiberMode.TWIST:
            active.append("ω (iso)" if cfg.isotropic_twist else "ω₁,ω₂,ω₃")
        elif cfg.fiber_mode == FiberMode.MIXING:  active.append("δ₀,δ₁,δ₂")
        elif cfg.fiber_mode == FiberMode.BOTH:    active.append("δ₀,δ₁,δ₂+ω₁,ω₂,ω₃")
        if getattr(cfg, 'enable_offdiag_shear', False): active.append("q₃,q₄,q₅")
        dof_str = ", ".join(active) if active else "none"
        self.logger.info(
            f"Topology: {cfg.topology.value.upper()}×S¹  "
            f"[Active DOF: {dof_str}]  "
            f"[Torsion: {self.mode.value}]"
        )
        self.logger.success("Setup complete")

    # =========================================================================
    # Fiber structure constant helpers — used by S³ and Nil³ subclasses
    # =========================================================================

    def _add_s3_twist_C(
        self,
        C: MutableDenseNDimArray,
        params: Dict,
        epsilon: Any,
        s: Any,
        r: Any,
        L: Any,
    ) -> None:
        """
        Add S³ twist structure constants C³_{jk} to C.

        Derived from ê³ = L(dτ + ω₁σ⁰ + ω₂σ¹ + ω₃σ²):
            C³_{12} = +4L ω₁ (1+ε)^{1/3}(1+s)   / r²
            C³_{20} = +4L ω₂ (1+ε)^{1/3}/(1+s)  / r²
            C³_{01} = +4L ω₃ (1+ε)^{−2/3}        / r²    (independent of s)
        """
        omega1, omega2, omega3 = params['omega1'], params['omega2'], params['omega3']

        twist_12 = (1 + epsilon) ** Rational(1, 3) * (1 + s)
        twist_20 = (1 + epsilon) ** Rational(1, 3) / (1 + s)
        twist_01 = (1 + epsilon) ** Rational(-2, 3)  # independent of s

        C[3, 1, 2] = +4 * L * omega1 * twist_12 / r ** 2
        C[3, 2, 1] = -4 * L * omega1 * twist_12 / r ** 2
        C[3, 2, 0] = +4 * L * omega2 * twist_20 / r ** 2
        C[3, 0, 2] = -4 * L * omega2 * twist_20 / r ** 2
        C[3, 0, 1] = +4 * L * omega3 * twist_01 / r ** 2
        C[3, 1, 0] = -4 * L * omega3 * twist_01 / r ** 2

    def _add_s3_mixing_C(
        self,
        C: MutableDenseNDimArray,
        params: Dict,
        factor_0: Any,
        factor_1: Any,
        factor_2: Any,
        r: Any,
        n_delta: Any,
    ) -> None:
        """
        Add S³–S¹ mixing structure constants to C.

        Derived from the mixing frame ansatz (rotating ê² and ê³ by δ):
            C⁰_{13} = +4δ · factor_0 / (r·n)
            C¹_{30} = +4δ · factor_1 / (r·n)
            C³_{01} = −4δ · factor_2 / (r·n)
        """
        delta = params['delta']

        C_mix_0 = 4 * delta * factor_0 / (r * n_delta)
        C[0, 1, 3] = +C_mix_0
        C[0, 3, 1] = -C_mix_0

        C_mix_1 = 4 * delta * factor_1 / (r * n_delta)
        C[1, 0, 3] = -C_mix_1
        C[1, 3, 0] = +C_mix_1

        C_mix_ax = 4 * delta * factor_2 / (r * n_delta)
        C[3, 0, 1] = -C_mix_ax
        C[3, 1, 0] = +C_mix_ax

    def _apply_s3_mixing_rotation_to_C(
        self,
        C: MutableDenseNDimArray,
        params: Dict,
    ) -> None:
        """
        Apply the 3-component mixing rotation M(δ₀,δ₁,δ₂) to all structure constants.

        M = R₀₃(δ₀) × R₁₃(δ₁) × R₂₃(δ₂)

        Each R_{i3} is an SO(2) rotation in the (ê^i, ê³) plane:
            tan θᵢ = δᵢ,  cos θᵢ = 1/nᵢ,  sin θᵢ = δᵢ/nᵢ,  nᵢ = √(1+δᵢ²)

        Tensor transformation rule:
            C'^a_{bc} = M^a_d (M^{-1})^e_b (M^{-1})^f_c C^d_{ef}

        M is orthogonal (subgroup of SO(4)), so M^{-1} = M^T.

        Limiting cases:
            δ₀=δ₁=δ₂=0  →  M=I, C unchanged  (reduces to TWIST-only)
            ω=0          →  reduces to MIXING-only structure constants
            δ₁=δ₂=0     →  matches old scalar δ implementation ((0,3)-plane version for δ₀ = δ)
            δ₀=δ₁=0     →  exactly matches old scalar δ implementation (δ₂ = δ) ✓
        """
        d0 = params.get('delta0', S.Zero)
        d1 = params.get('delta1', S.Zero)
        d2 = params.get('delta2', S.Zero)

        # No rotation needed when all δ are Zero (fast path)
        if d0 == S.Zero and d1 == S.Zero and d2 == S.Zero:
            return

        # Use polynomial auxiliary symbols: cd_i = 1/√(1+δ_i²), sd_i = δ_i/√(1+δ_i²)
        # Use cd/sd symbols from params if available (pipeline speedup).
        # Fall back to direct sqrt computation otherwise.
        def _cs(d, cd_key, sd_key):
            cd = params.get(cd_key)
            sd = params.get(sd_key)
            if cd is not None and sd is not None:
                # When δ is Zero, use identity transformation values directly
                if d == S.Zero:
                    return S.One, S.Zero
                return cd, sd
            # Fallback: direct sqrt computation (legacy behavior)
            n = sqrt(1 + d ** 2)
            return 1 / n, d / n

        c0, s0 = _cs(d0, 'cd0', 'sd0')
        c1, s1 = _cs(d1, 'cd1', 'sd1')
        c2, s2 = _cs(d2, 'cd2', 'sd2')

        # ── Build 4×4 rotation matrix ────────────────────────────────────────────
        # R₀₃: rotation in the (ê⁰, ê³) plane
        #   ê'⁰ = c0 ê⁰ + s0 ê³
        #   ê'³ = −s0 ê⁰ + c0 ê³
        R03 = Matrix([
            [ c0, S.Zero, S.Zero, s0],
            [S.Zero, S.One, S.Zero, S.Zero],
            [S.Zero, S.Zero, S.One, S.Zero],
            [-s0, S.Zero, S.Zero, c0],
        ])

        # R₁₃: rotation in the (ê¹, ê³) plane
        R13 = Matrix([
            [S.One, S.Zero, S.Zero, S.Zero],
            [S.Zero,  c1, S.Zero, s1],
            [S.Zero, S.Zero, S.One, S.Zero],
            [S.Zero, -s1, S.Zero, c1],
        ])

        # R₂₃: rotation in the (ê², ê³) plane  ← equivalent to old scalar δ implementation
        R23 = Matrix([
            [S.One, S.Zero, S.Zero, S.Zero],
            [S.Zero, S.One, S.Zero, S.Zero],
            [S.Zero, S.Zero,  c2, s2],
            [S.Zero, S.Zero, -s2, c2],
        ])

        M  = R03 * R13 * R23   # composite rotation (orthogonal matrix)
        # In the tensor transformation, Mi[e,b] = (M^{-1})[b,e] = M^T[b,e] = M[e,b].
        # Therefore for an orthogonal matrix, Mi = M is directly the correct index ordering.
        Mi = M

        dim = 4

        # ── Collect nonzero old entries ─────────────────────────────────────────
        C_old = {}
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    v = C[a, b, c]
                    if v != S.Zero:
                        C_old[(a, b, c)] = v

        # ── Zero-clear C ───────────────────────────────────────────────
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    C[a, b, c] = S.Zero

        # ── Apply tensor transformation ─────────────────────────────────────────
        # C'[a,b,c] = Σ_{d,e,f} M[a,d] × Mi[e,b] × Mi[f,c] × C_old[d,e,f]
        # Mi[e,b] = M[e,b] represents (M^{-1})[b,e] (index ordering from orthogonal Mi = M)
        for (d, e, f), v_old in C_old.items():
            for a in range(dim):
                ma_d = M[a, d]
                if ma_d == S.Zero:
                    continue
                for b in range(dim):
                    mi_e_b = Mi[e, b]
                    if mi_e_b == S.Zero:
                        continue
                    for c in range(dim):
                        mi_f_c = Mi[f, c]
                        if mi_f_c == S.Zero:
                            continue
                        C[a, b, c] = C[a, b, c] + ma_d * mi_e_b * mi_f_c * v_old

    # =========================================================================
    # Analysis helpers — public API called from analysis scripts
    # =========================================================================

    def get_riemann_lambdified(self, source=None):
        """
        Lambdify R_{abcd} and return a function usable for numerical Pontryagin scanning.
        If source is omitted, the engine's configured riemann_source is used.

        Returns
        -------
        func : callable
            f(*args) → numpy array (dim, dim, dim, dim)
        arg_names : list[str]
            ordered list of argument names for func
        """
        from sympy import lambdify, S as _S
        import numpy as np

        source_key = self._normalize_curvature_source(
            source if source is not None else self.riemann_source
        )
        tensor_key = 'riemann_abcd_LC' if source_key == 'lc' else 'riemann_abcd'

        if tensor_key not in self.data:
            raise RuntimeError(
                "Run engine first (call engine.run() before get_riemann_lambdified())."
            )

        params = self.data['params']
        R_sym  = self.data[tensor_key]
        dim    = self.data['dim']

        # Candidate list covering both S³ 'r' and Nil³ 'R'
        # Placing 'R' first ensures it is selected for Nil³ (where R is an alias for r)
        candidates = [
            'R', 'r', 'omega', 'omega1', 'omega2', 'omega3',
            'eta', 'V', 'L', 'kappa', 'theta_NY', 'alpha',
            'epsilon', 's', 'sigma',
            'delta0', 'delta1', 'delta2',
            'cd0', 'sd0', 'cd1', 'sd1', 'cd2', 'sd2',
            'q3', 'q4', 'q5',
        ]
        arg_symbols, arg_names = [], []
        seen_syms: set = set()
        for name in candidates:
            sym = params.get(name, _S.Zero)
            if hasattr(sym, 'is_Symbol') and sym.is_Symbol and sym not in seen_syms:
                seen_syms.add(sym)
                arg_symbols.append(sym)
                arg_names.append(name)

        component_funcs = {}
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):
                        expr = R_sym[a, b, c, d]
                        if expr != _S.Zero:
                            component_funcs[(a, b, c, d)] = lambdify(
                                arg_symbols, expr, modules='numpy'
                            )

        def R_numerical(*args):
            R = np.zeros((dim, dim, dim, dim))
            for (a, b, c, d), f in component_funcs.items():
                R[a, b, c, d] = float(f(*args))
            return R

        return R_numerical, arg_names

    def get_free_params(self) -> Dict[str, Any]:
        """
        Return a dict of currently active (non-Zero) SymPy Symbols.

        Returns
        -------
        dict[str, Symbol]
        """
        return {
            k: v
            for k, v in self.data['params'].items()
            if hasattr(v, 'is_Symbol') and v.is_Symbol
        }

    def get_mixing_cs_subs(self) -> Dict:
        """
        Return a substitution dict mapping polynomial auxiliary symbols (cd_i, sd_i)
        to their actual expressions in δ.

        The pipeline treats cd_i = cos θ_i, sd_i = sin θ_i as formal polynomial
        symbols for performance. Apply this substitution dict before numerical evaluation.

        Returns
        -------
        dict : {cd0: 1/sqrt(1+delta0**2), sd0: delta0/sqrt(1+delta0**2), ...}
               Entries for δ_i = Zero are omitted.

        Usage
        -----
        subs = engine.get_mixing_cs_subs()
        Veff_phys = engine.data['potential'].subs(subs)
        # Then substitute numerical values such as delta0=0, delta1=0, delta2=val
        """
        params = self.data['params']
        result = {}
        for i, name in enumerate(['0', '1', '2']):
            d = params.get(f'delta{name}', S.Zero)
            cd = params.get(f'cd{name}')
            sd = params.get(f'sd{name}')
            if cd is not None and hasattr(cd, 'is_Symbol') and cd.is_Symbol:
                if d == S.Zero:
                    result[cd] = S.One
                    result[sd] = S.Zero
                else:
                    n = sqrt(1 + d ** 2)
                    result[cd] = 1 / n
                    result[sd] = d / n
        return result

    def get_mixing_cs_subs_numerical(self, d0_val: float = 0.0, d1_val: float = 0.0, d2_val: float = 0.0) -> Dict:
        """
        Return a dict substituting polynomial auxiliary symbols (cd_i, sd_i) with numerical values.

        Parameters
        ----------
        d0_val, d1_val, d2_val : float
            Numerical values of δ₀, δ₁, δ₂.

        Returns
        -------
        dict : {cd0: float, sd0: float, cd1: float, ...}
        """
        import math
        params = self.data['params']
        result = {}
        for name, d_val in [('0', d0_val), ('1', d1_val), ('2', d2_val)]:
            cd = params.get(f'cd{name}')
            sd = params.get(f'sd{name}')
            if cd is not None and hasattr(cd, 'is_Symbol') and cd.is_Symbol:
                n = math.sqrt(1.0 + d_val ** 2)
                result[cd] = 1.0 / n
                result[sd] = d_val / n
        return result

    def get_config_summary(self) -> str:
        """Return a one-line summary of the configuration."""
        from .unified import FiberMode
        cfg   = self.config
        parts = [cfg.topology.value.upper()]
        if cfg.enable_squash:  parts.append("squash(ε)")
        if cfg.enable_shear:   parts.append("shear(s)")
        if cfg.fiber_mode == FiberMode.TWIST:
            parts.append("twist(ω)" if cfg.isotropic_twist else "twist(ω₁ω₂ω₃)")
        elif cfg.fiber_mode == FiberMode.MIXING:
            parts.append("mixing(δ₀δ₁δ₂)")
        elif cfg.fiber_mode == FiberMode.BOTH:
            parts.append("both(δ₀δ₁δ₂,ω₁ω₂ω₃)")
        if getattr(cfg, 'enable_offdiag_shear', False):
            parts.append("offdiag_shear(q₃q₄q₅)")
        parts.append(f"Weyl:{cfg.weyl_source.value.upper()}")
        parts.append(f"Riemann:{cfg.riemann_source.value.upper()}")
        parts.append(self.mode.value)
        return " + ".join(parts)
