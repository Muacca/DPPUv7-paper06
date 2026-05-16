"""
Scan SD Diagnostics
===================

Self-Duality diagnostics integrated with scan results.
Evaluates SD conditions at r = r*(V, eta) from parameter scan.
"""

import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import minimize_scalar

from .scan_results_loader import ScanResultsLoader
from ..curvature import CurvatureSDDiagnostics


@dataclass
class SDScanResult:
    """Result container for SD parameter scan."""
    eta: float
    V: float
    r_star: float
    sd_residual: float
    asd_residual: float
    curvature_norm: float
    is_sd: bool
    is_asd: bool
    is_nontrivial: bool
    is_nontrivial_sd: bool
    phase_type: str


class SDScanDiagnostics:
    """
    Self-Duality diagnostics integrated with scan results.

    Key features:
    1. Evaluates SD conditions at r = r*(V, eta) from scan results
    2. Cross-references with stability type classification
    3. Identifies SD curves in stable (Type I/II) regions

    Usage:
        from dppu.topology import CurvatureSource, TopologyType, make_engine
        from dppu.torsion import Mode, NyVariant
        from dppu.curvature import SDExtensionMixin
        from dppu.scanning import ScanResultsLoader, SDScanDiagnostics

        engine = make_engine(
            TopologyType.S3, Mode.MX, NyVariant.FULL, enable_squash=True,
            weyl_source=CurvatureSource.LC,
            riemann_source=CurvatureSource.LC,
        )
        engine.run()
        SDExtensionMixin.attach_to(engine)

        loader = ScanResultsLoader.from_csv('scan.csv', theta_fixed=0.0)
        sd_diag = SDScanDiagnostics(engine, loader)
        results = sd_diag.scan_with_phase1_rstar(...)
    """

    def __init__(
        self,
        engine,
        scan_loader: Optional[ScanResultsLoader] = None
    ):
        """
        Initialize with engine and optional scan results.

        Args:
            engine: BaseFrameEngine with SDExtensionMixin attached
            scan_loader: ScanResultsLoader instance
        """
        self.engine = engine
        self.base_diag = CurvatureSDDiagnostics(engine)
        self.scan_loader = scan_loader

    def get_r_star(self, V: float, eta: float, theta_NY: float = 0.0) -> Optional[float]:
        """Get stable radius r* from scan results."""
        if self.scan_loader is not None:
            return self.scan_loader.get_r_star(V, eta)

        warnings.warn("No scan results loaded. Using fallback r=1.0")
        return 1.0

    def get_phase_type(self, V: float, eta: float) -> str:
        """Get stability classification (Type I/II/III)."""
        if self.scan_loader is not None:
            return self.scan_loader.get_phase_type(V, eta)
        return "unknown"

    def evaluate_at_rstar(
        self,
        V: float,
        eta: float,
        L: float = 1.0,
        kappa: float = 1.0,
        theta_NY: float = 0.0,
        eps_sd: float = 1e-6,
        eps_R: float = 1e-8,
        extra_params: Optional[Dict] = None
    ) -> Optional[SDScanResult]:
        """
        Evaluate SD status at the Phase 1 stable point r*.

        Args:
            V, eta: Torsion parameters
            L, kappa, theta_NY: Other parameters
            eps_sd, eps_R: Thresholds
            extra_params: Additional topology-specific parameters

        Returns:
            SDScanResult or None if no stable point exists
        """
        r_star = self.get_r_star(V, eta, theta_NY)
        if r_star is None:
            return None

        params = {
            'r': r_star,
            'L': L,
            'eta': eta,
            'V': V,
            'kappa': kappa,
            'theta_NY': theta_NY
        }

        if extra_params:
            params.update(extra_params)

        result = self.base_diag.evaluate_sd_status(params, eps_sd, eps_R)
        phase_type = self.get_phase_type(V, eta)

        return SDScanResult(
            eta=eta,
            V=V,
            r_star=r_star,
            sd_residual=result['sd_residual'],
            asd_residual=result['asd_residual'],
            curvature_norm=result['curvature_norm'],
            is_sd=result['is_sd'],
            is_asd=result['is_asd'],
            is_nontrivial=result['is_nontrivial'],
            is_nontrivial_sd=result['is_nontrivial_sd'],
            phase_type=phase_type
        )

    def scan_with_phase1_rstar(
        self,
        eta_range: Tuple[float, float, int],
        V_range: Tuple[float, float, int],
        L: float = 1.0,
        kappa: float = 1.0,
        theta_NY: float = 0.0,
        eps_sd: float = 1e-6,
        eps_R: float = 1e-8,
        verbose: bool = True,
        extra_params: Optional[Dict] = None
    ) -> Dict:
        """
        Scan (eta, V) plane using r* from Phase 1 at each point.

        Args:
            eta_range, V_range: (min, max, n_points)
            L, kappa, theta_NY: Fixed parameters
            eps_sd, eps_R: Thresholds
            verbose: Print progress
            extra_params: Additional topology-specific parameters

        Returns:
            Dict with results, sd_curve, asd_curve, type intersections
        """
        eta_arr = np.linspace(eta_range[0], eta_range[1], eta_range[2])
        V_arr = np.linspace(V_range[0], V_range[1], V_range[2])

        results = []
        sd_curve = []
        asd_curve = []
        type_I_sd = []
        type_II_sd = []

        total = len(eta_arr) * len(V_arr)
        count = 0

        for eta in eta_arr:
            for V in V_arr:
                count += 1
                if verbose and count % 100 == 0:
                    print(f"Progress: {count}/{total} ({100*count/total:.1f}%)")

                result = self.evaluate_at_rstar(
                    V=V, eta=eta, L=L, kappa=kappa, theta_NY=theta_NY,
                    eps_sd=eps_sd, eps_R=eps_R, extra_params=extra_params
                )

                if result is None:
                    continue

                results.append(result)

                if result.is_nontrivial_sd:
                    if result.is_sd:
                        sd_curve.append((eta, V))
                        if result.phase_type == "I":
                            type_I_sd.append((eta, V))
                        elif result.phase_type == "II":
                            type_II_sd.append((eta, V))
                    if result.is_asd:
                        asd_curve.append((eta, V))

        return {
            'results': results,
            'sd_curve': sd_curve,
            'asd_curve': asd_curve,
            'type_I_sd_intersection': type_I_sd,
            'type_II_sd_intersection': type_II_sd,
            'scan_params': {
                'eta_range': eta_range,
                'V_range': V_range,
                'L': L,
                'kappa': kappa,
                'theta_NY': theta_NY,
                'eps_sd': eps_sd,
                'eps_R': eps_R
            }
        }

    def find_sd_curve_by_minimization(
        self,
        eta_range: Tuple[float, float],
        V_fixed: float,
        L: float = 1.0,
        kappa: float = 1.0,
        theta_NY: float = 0.0,
        n_initial: int = 20
    ) -> List[Tuple[float, float, str]]:
        """
        Find SD curve by minimizing min(sd_res, asd_res) along eta.

        Args:
            eta_range: (min, max)
            V_fixed: Fixed V value for this slice

        Returns:
            List of (eta, residual, type) where type is 'SD' or 'ASD'
        """
        def objective(eta):
            result = self.evaluate_at_rstar(
                V=V_fixed, eta=eta, L=L, kappa=kappa, theta_NY=theta_NY
            )
            if result is None:
                return 1e10
            return min(result.sd_residual, result.asd_residual)

        eta_arr = np.linspace(eta_range[0], eta_range[1], n_initial)
        candidates = []

        for i in range(len(eta_arr) - 1):
            res1 = objective(eta_arr[i])
            res2 = objective(eta_arr[i+1])

            if res1 < 1e-3 or res2 < 1e-3:
                result = minimize_scalar(
                    objective,
                    bounds=(eta_arr[i], eta_arr[i+1]),
                    method='bounded'
                )

                if result.fun < 1e-6:
                    final = self.evaluate_at_rstar(
                        V=V_fixed, eta=result.x, L=L, kappa=kappa, theta_NY=theta_NY
                    )
                    if final:
                        sd_type = 'SD' if final.sd_residual < final.asd_residual else 'ASD'
                        candidates.append((result.x, result.fun, sd_type))

        return candidates
