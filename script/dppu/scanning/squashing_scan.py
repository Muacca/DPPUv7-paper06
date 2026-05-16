"""
Squashing Parameter Scan
========================

Scans the effective potential V_eff(r, epsilon) to visualize the landscape.
"""

import numpy as np
from typing import Dict, Any, Tuple, Optional

from ..torsion.mode import Mode
from ..torsion.nieh_yan import NyVariant
from ..topology import CurvatureSource, TopologyType, make_engine


def get_engine(topology: str, mode: Mode, ny_variant: NyVariant):
    """Factory for topology engines."""
    if topology == 'S3':
        return make_engine(
            TopologyType.S3, mode, ny_variant, enable_squash=True,
            weyl_source=CurvatureSource.LC,
            riemann_source=CurvatureSource.LC,
        )
    elif topology == 'T3':
        return make_engine(
            TopologyType.T3, mode, ny_variant,
            weyl_source=CurvatureSource.LC,
            riemann_source=CurvatureSource.LC,
        )
    elif topology == 'Nil3':
        return make_engine(
            TopologyType.NIL3, mode, ny_variant, enable_squash=True,
            weyl_source=CurvatureSource.LC,
            riemann_source=CurvatureSource.LC,
        )
    else:
        raise ValueError(f"Unknown topology: {topology}")


def scan_potential_landscape(
    topology: str,
    mode: str,  # "MX", "AX", or "VT"
    ny_variant: str, # "FULL", "TT", "REE"
    params: Dict[str, float],
    r_range: np.ndarray,
    epsilon_range: np.ndarray
) -> Dict[str, Any]:
    """
    Compute V_eff on a (r, epsilon) grid.

    Args:
        topology: "S3", "T3", or "Nil3"
        mode: Mode string
        ny_variant: NyVariant string
        params: Dictionary containing fixed parameters:
                {'V', 'eta', 'theta_NY', 'L', 'kappa', 'alpha'}
        r_range: 1D array of r values
        epsilon_range: 1D array of epsilon values

    Returns:
        Dictionary with keys: 'R', 'Epsilon', 'V_eff', 'params'
    """
    # Convert strings to enums
    mode_enum = Mode[mode]
    ny_enum = NyVariant[ny_variant]
    
    # Initialize and run engine
    engine = get_engine(topology, mode_enum, ny_enum)
    # We need to run the engine to define symbols and potential
    # The run() method executes the pipeline purely symbolically
    engine.run()
    
    # Get numerical function
    # V_func(r, V, eta, theta, L, kappa, epsilon, alpha)
    V_func = engine.get_effective_potential_function()
    
    # Prepare grid
    R, Eps = np.meshgrid(r_range, epsilon_range)
    
    # Extract parameters
    # Default values compatible with scripts
    V_val = params.get('V', 1.0)
    eta_val = params.get('eta', 0.0)
    theta_val = params.get('theta_NY', 0.0)
    L_val = params.get('L', 1.0)
    kappa_val = params.get('kappa', 1.0)
    alpha_val = params.get('alpha', 0.0)
    
    # Evaluate potential
    # Note: engine.pipeline.py get_effective_potential_function args order:
    # r, V, eta, theta_NY, L, kappa, epsilon, alpha
    
    # We pass meshgrids R and Eps, other args are scalars.
    # numpy broadcasting should handle this.
    Z = V_func(R, V_val, eta_val, theta_val, L_val, kappa_val, Eps, alpha_val)
    
    # Replace large values or NaNs if necessary?
    # Potentials can be singular at r=0. V_func might return inf.
    # We accept what numpy returns.
    
    return {
        'R': R,
        'Epsilon': Eps,
        'V_eff': Z,
        'params': params,
        'topology': topology,
        'mode': mode,
        'ny_variant': ny_variant
    }
