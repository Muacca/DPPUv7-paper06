"""
Self-Duality Extension
======================

Provides R^{ab}_{cd} computation and SD diagnostics integration.
"""

import types
from typing import Any, Callable, Dict, Optional, Tuple

import numpy as np
from sympy import S, cancel, lambdify
from sympy.tensor.array import MutableDenseNDimArray


class SDExtensionMixin:
    """
    Mixin to add Self-Duality methods to BaseFrameEngine.

    Usage:
        engine = make_engine(TopologyType.S3, mode, variant, enable_squash=True)
        engine.run()
        SDExtensionMixin.attach_to(engine)
        R = engine.get_R_ab_cd_numerical(params)
    """

    _R_ab_cd_lambdified: Optional[Callable] = None
    _param_symbols: Optional[Dict] = None

    @classmethod
    def attach_to(cls, engine, source=None):
        """Dynamically attach SD methods to engine instance."""
        source_key = engine._normalize_curvature_source(
            source if source is not None else getattr(engine, "riemann_source", "ec")
        )
        engine.get_R_ab_cd_numerical = types.MethodType(
            cls.get_R_ab_cd_numerical, engine
        )
        engine._compile_R_ab_cd_lambdified = types.MethodType(
            cls._compile_R_ab_cd_lambdified, engine
        )
        engine._build_R_ab_cd_symbolic = types.MethodType(
            cls._build_R_ab_cd_symbolic, engine
        )
        engine._get_unified_param_symbols = types.MethodType(
            cls._get_unified_param_symbols, engine
        )

        engine._R_ab_cd_lambdified = None
        engine._R_ab_cd_symbolic = None
        engine._param_symbols_unified = None
        engine._sd_curvature_source = source_key

        return engine

    def _get_unified_param_symbols(self) -> Tuple[Dict, list]:
        """Get unified parameter symbols."""
        if hasattr(self, '_param_symbols_unified') and self._param_symbols_unified:
            return self._param_symbols_unified

        params = self.data.get('params', {})
        unified = {}

        for key in ['r', 'L', 'kappa', 'theta_NY']:
            if key in params:
                unified[key] = params[key]
            elif key == 'r' and 'R1' in params:
                unified['r'] = params['R1']
            elif key == 'r' and 'R' in params:
                unified['r'] = params['R']

        for key in ['eta', 'V']:
            unified[key] = params.get(key, S.Zero)

        symbol_list = []
        for name in ['r', 'L', 'eta', 'V', 'kappa', 'theta_NY']:
            sym = unified.get(name, S.Zero)
            if hasattr(sym, 'is_Symbol') and sym.is_Symbol:
                symbol_list.append(sym)

        # Include any additional symbolic params (e.g. epsilon, alpha)
        seen = {str(s) for s in symbol_list}
        for name, sym in params.items():
            if name in unified:
                continue
            if hasattr(sym, 'is_Symbol') and sym.is_Symbol and str(sym) not in seen:
                unified[str(sym)] = sym
                symbol_list.append(sym)
                seen.add(str(sym))
            elif hasattr(sym, 'free_symbols') and sym.free_symbols:
                for s in sym.free_symbols:
                    if str(s) not in seen:
                        unified[str(s)] = s
                        symbol_list.append(s)
                        seen.add(str(s))

        self._param_symbols_unified = (unified, symbol_list)
        return unified, symbol_list

    def _build_R_ab_cd_symbolic(self) -> MutableDenseNDimArray:
        """Build R^{ab}_{cd} from R^a_{bcd}."""
        if hasattr(self, '_R_ab_cd_symbolic') and self._R_ab_cd_symbolic:
            return self._R_ab_cd_symbolic

        source_key = getattr(self, "_sd_curvature_source", "ec")
        riemann_key = 'riemann_LC' if source_key == 'lc' else 'riemann'
        R_a_bcd = self.data[riemann_key]
        eta = self.data['metric_frame']
        dim = self.data['dim']

        R_ab_cd = MutableDenseNDimArray.zeros(dim, dim, dim, dim)

        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):
                        val = S.Zero
                        for e in range(dim):
                            val += eta[b, e] * R_a_bcd[a, e, c, d]
                        if val != S.Zero:
                            R_ab_cd[a, b, c, d] = cancel(val)

        self._R_ab_cd_symbolic = R_ab_cd
        return R_ab_cd

    def _compile_R_ab_cd_lambdified(self) -> Callable:
        """Compile R^{ab}_{cd} for fast numerical evaluation."""
        if hasattr(self, '_R_ab_cd_lambdified') and self._R_ab_cd_lambdified:
            return self._R_ab_cd_lambdified

        R_ab_cd = self._build_R_ab_cd_symbolic()
        _, symbol_list = self._get_unified_param_symbols()
        dim = self.data['dim']

        component_funcs = {}
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):
                        expr = R_ab_cd[a, b, c, d]
                        if expr != S.Zero:
                            if symbol_list:
                                component_funcs[(a, b, c, d)] = lambdify(
                                    symbol_list, expr, modules='numpy'
                                )
                            else:
                                val = float(expr)
                                component_funcs[(a, b, c, d)] = lambda *args, v=val: v

        def R_numerical(params_dict: Dict[str, float]) -> np.ndarray:
            result = np.zeros((dim, dim, dim, dim))
            args = [params_dict.get(str(sym), 0.0) for sym in symbol_list]
            for (a, b, c, d), func in component_funcs.items():
                result[a, b, c, d] = float(func(*args))
            return result

        self._R_ab_cd_lambdified = R_numerical
        return R_numerical

    def get_R_ab_cd_numerical(self, params_dict: Dict[str, float]) -> np.ndarray:
        """Get curvature tensor R^{ab}_{cd} as numpy array."""
        R_func = self._compile_R_ab_cd_lambdified()
        return R_func(params_dict)
