"""
Pontryagin Diagnostics
======================

SD/ASD analysis and Pontryagin inner product P = <R, *R>.
"""

from typing import Dict
import numpy as np
from .hodge import compute_hodge_dual


class CurvatureSDDiagnostics:
    """
    Self-Duality Diagnostics for Curvature Tensor.

    Computes:
    - E = <R, R>
    - P = <R, *R> (Pontryagin inner product)
    - SD/ASD residuals
    """

    def __init__(self, engine):
        """Initialize with engine that has SD extension attached."""
        self.engine = engine
        if not hasattr(engine, 'get_R_ab_cd_numerical'):
            raise RuntimeError("Engine needs SDExtensionMixin attached.")

    def compute_hodge_dual(self, R: np.ndarray) -> np.ndarray:
        """Delegate to module-level compute_hodge_dual."""
        return compute_hodge_dual(R)

    def evaluate_sd_status(
        self,
        params_dict: Dict[str, float],
        eps_sd: float = 1e-6,
        eps_R: float = 1e-8
    ) -> Dict:
        """Evaluate SD/ASD status at parameter point."""
        tiny = 1e-30

        R = self.engine.get_R_ab_cd_numerical(params_dict)
        R_dual = compute_hodge_dual(R)

        R_flat = R.ravel()
        R_dual_flat = R_dual.ravel()

        E_RR = float(np.dot(R_flat, R_flat))
        P_RstarR = float(np.dot(R_flat, R_dual_flat))
        p_ratio = P_RstarR / max(E_RR, tiny)

        Rplus_norm2 = 0.5 * (E_RR + P_RstarR)
        Rminus_norm2 = 0.5 * (E_RR - P_RstarR)
        Rplus_frac = Rplus_norm2 / max(E_RR, tiny)
        Rminus_frac = Rminus_norm2 / max(E_RR, tiny)

        diff_sd = R_flat - R_dual_flat
        sd_residual2_numeric = float(np.dot(diff_sd, diff_sd))
        sd_residual = np.sqrt(sd_residual2_numeric)
        sd_residual2_formula = 2.0 * (E_RR - P_RstarR)
        sd_residual2_delta = abs(sd_residual2_numeric - sd_residual2_formula)

        sum_asd = R_flat + R_dual_flat
        asd_residual2_numeric = float(np.dot(sum_asd, sum_asd))
        asd_residual = np.sqrt(asd_residual2_numeric)
        asd_residual2_formula = 2.0 * (E_RR + P_RstarR)
        asd_residual2_delta = abs(asd_residual2_numeric - asd_residual2_formula)

        curvature_norm = np.sqrt(E_RR)

        is_sd = sd_residual < eps_sd
        is_asd = asd_residual < eps_sd
        is_nontrivial = curvature_norm > eps_R

        return {
            'sd_residual': sd_residual,
            'asd_residual': asd_residual,
            'curvature_norm': curvature_norm,
            'is_sd': is_sd,
            'is_asd': is_asd,
            'is_nontrivial': is_nontrivial,
            'is_nontrivial_sd': (is_sd or is_asd) and is_nontrivial,
            'params': params_dict.copy(),
            'E_RR': E_RR,
            'P_RstarR': P_RstarR,
            'p_P_over_E': p_ratio,
            'Rplus_norm2': Rplus_norm2,
            'Rminus_norm2': Rminus_norm2,
            'Rplus_frac': Rplus_frac,
            'Rminus_frac': Rminus_frac,
            'sd_residual2_numeric': sd_residual2_numeric,
            'sd_residual2_formula': sd_residual2_formula,
            'sd_residual2_delta': sd_residual2_delta,
            'asd_residual2_numeric': asd_residual2_numeric,
            'asd_residual2_formula': asd_residual2_formula,
            'asd_residual2_delta': asd_residual2_delta,
        }

    def scan_parameter_plane(
        self,
        r_val: float,
        L_val: float,
        kappa_val: float,
        theta_NY_val: float,
        eta_range: tuple,
        V_range: tuple,
        eps_sd: float = 1e-6,
        eps_R: float = 1e-8
    ) -> Dict:
        """Scan (eta, V) plane for SD conditions."""
        eta_arr = np.linspace(eta_range[0], eta_range[1], eta_range[2])
        V_arr = np.linspace(V_range[0], V_range[1], V_range[2])

        sd_map = np.zeros((len(eta_arr), len(V_arr)))
        sd_points = []

        for i, eta in enumerate(eta_arr):
            for j, V in enumerate(V_arr):
                params = {
                    'r': r_val, 'L': L_val, 'eta': eta, 'V': V,
                    'kappa': kappa_val, 'theta_NY': theta_NY_val
                }
                result = self.evaluate_sd_status(params, eps_sd, eps_R)
                sd_map[i, j] = result['sd_residual']
                if result['is_nontrivial_sd']:
                    sd_points.append((eta, V))

        return {
            'eta_grid': eta_arr,
            'V_grid': V_arr,
            'sd_residual_map': sd_map,
            'sd_points': sd_points,
        }
