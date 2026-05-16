"""
EC Action Builder
=================

Helper functions for building Einstein-Cartan + Weyl-squared effective potential.

These consolidate the 6-line C²_EC boilerplate that appeared in every phase script.

Functions:
- compute_c2_ec: C²_EC from engine data dict
- build_veff_ec: symbolic V_eff_EC and V_eff_LC (with caching)
- build_veff_ec_func: lambdified V_eff_EC(r, eta, V, alpha, theta)
"""

from typing import Any, Callable, Dict, Optional, Tuple

from sympy import S, cancel, lambdify

from ..curvature.ricci import compute_ricci_tensor
from ..curvature.weyl import compute_weyl_tensor, compute_weyl_scalar
from .lagrangian import compute_lagrangian, compute_action
from .potential import compute_effective_potential

# --------------------------------------------------------------------------
# Internal cache keyed by (topology_type, torsion_mode, ny_variant, extras)
# --------------------------------------------------------------------------
_veff_cache: Dict[tuple, Any] = {}


def compute_c2_ec(data: Dict) -> Any:
    """
    Compute C²_EC (Weyl scalar for EC connection) from engine.data.

    Consolidates the repeated 3-step pattern:
        ricci_EC = compute_ricci_tensor(data['riemann'], data['dim'])
        weyl_EC  = compute_weyl_tensor(data['riemann_abcd'], ricci_EC, ...)
        c2_ec    = cancel(compute_weyl_scalar(weyl_EC, metric_inv, data['dim']))

    Args:
        data: dict returned by UnifiedEngine after engine.run()

    Returns:
        SymPy expression C²_EC (symbolic)
    """
    ricci_EC = compute_ricci_tensor(data['riemann'], data['dim'])
    metric    = data['metric_frame']
    weyl_EC   = compute_weyl_tensor(
        data['riemann_abcd'], ricci_EC,
        data['ricci_scalar'], metric, data['dim']
    )
    return cancel(compute_weyl_scalar(weyl_EC, metric.inv(), data['dim']))


def build_veff_ec(
    topology_type: Any,
    torsion_mode: Any = None,
    ny_variant: Any = None,
    fiber_mode: Any = None,
    enable_squash: bool = False,
    enable_shear: bool = False,
    L_val: Optional[float] = None,
    kappa_val: Optional[float] = None,
    use_cache: bool = True,
) -> Tuple[Any, Any, Dict]:
    """
    Build V_eff_EC and V_eff_LC symbolically using UnifiedEngine.

    Args:
        topology_type: TopologyType enum value (S3, T3, NIL3)
        torsion_mode: Mode enum (default Mode.MX)
        ny_variant: NyVariant enum (default NyVariant.FULL)
        fiber_mode: FiberMode enum (default FiberMode.NONE)
        enable_squash: enable squash DOF (default False)
        enable_shear: enable shear DOF (default False)
        L_val: if given, substitute L = L_val (default: leave symbolic)
        kappa_val: if given, substitute kappa = kappa_val (default: leave symbolic)
        use_cache: use module-level cache (default True)

    Returns:
        (veff_ec, veff_lc, params) — all symbolic (with optional L/κ substituted)
    """
    # Inline imports avoid circular dependency at module load time
    from ..topology.unified import UnifiedEngine, DOFConfig, FiberMode as FM
    from ..torsion.mode import Mode
    from ..torsion.nieh_yan import NyVariant as NyV

    if torsion_mode is None:
        torsion_mode = Mode.MX
    if ny_variant is None:
        ny_variant = NyV.FULL
    if fiber_mode is None:
        fiber_mode = FM.NONE

    cache_key = (topology_type, torsion_mode, ny_variant, fiber_mode,
                 enable_squash, enable_shear, L_val, kappa_val)
    if use_cache and cache_key in _veff_cache:
        return _veff_cache[cache_key]

    cfg = DOFConfig(
        topology=topology_type,
        torsion_mode=torsion_mode,
        fiber_mode=fiber_mode,
        enable_squash=enable_squash,
        enable_shear=enable_shear,
        ny_variant=ny_variant,
    )
    eng = UnifiedEngine(cfg)
    eng.run()
    data = eng.data
    params = data['params']

    # V_eff_EC: pipeline now stores C²_EC in data['potential'] (step E4.7b → E4.11)
    veff_ec = data['potential']

    # V_eff_LC: recompute from C²_LC explicitly so callers can compare LC vs EC.
    # data['weyl_scalar'] = C²_LC (stored by pipeline step E4.3b, unchanged).
    c2_lc = data.get('weyl_scalar', S.Zero)
    L_lc = compute_lagrangian(
        data['ricci_scalar'], data['nieh_yan_density'],
        params['kappa'], params['theta_NY'],
        weyl_scalar=c2_lc, alpha=params['alpha'],
    )
    veff_lc = compute_effective_potential(compute_action(L_lc, data['total_volume']))

    if L_val is not None or kappa_val is not None:
        subs = {}
        if L_val is not None:
            subs[params['L']] = L_val
        if kappa_val is not None:
            subs[params['kappa']] = kappa_val
        veff_ec = veff_ec.subs(subs)
        veff_lc = veff_lc.subs(subs)

    result = (veff_ec, veff_lc, params)
    if use_cache:
        _veff_cache[cache_key] = result
    return result


def build_veff_ec_func(
    topology_type: Any,
    torsion_mode: Any = None,
    ny_variant: Any = None,
    fiber_mode: Any = None,
    enable_squash: bool = False,
    enable_shear: bool = False,
) -> Callable:
    """
    Build a lambdified V_eff_EC(r, eta, V, alpha, theta) with L=kappa=1.

    Args:
        topology_type, torsion_mode, ny_variant, fiber_mode,
        enable_squash, enable_shear: passed to build_veff_ec

    Returns:
        Callable f(r, eta, V, alpha, theta) → float
    """
    veff_ec, _, params = build_veff_ec(
        topology_type,
        torsion_mode=torsion_mode,
        ny_variant=ny_variant,
        fiber_mode=fiber_mode,
        enable_squash=enable_squash,
        enable_shear=enable_shear,
        L_val=1.0,
        kappa_val=1.0,
    )

    syms = [params['r'], params['eta'], params['V'],
            params['alpha'], params['theta_NY']]
    return lambdify(syms, veff_ec, modules='numpy')
