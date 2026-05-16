"""
Scan Results Loader
===================

Loader and interpolator for parameter scan results.
"""

import warnings
from pathlib import Path
from typing import Dict, Optional, Union

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator


class ScanResultsLoader:
    """
    Loader and interpolator for parameter scan results.

    Features:
    - Loads CSV files from parameter_scan
    - 2D interpolation on (eta, V) plane for fixed theta
    - Handles missing values (type-III points)

    Usage:
        loader = ScanResultsLoader.from_csv(
            'output/dppu_scan_S3_FULL.csv',
            theta_fixed=0.0
        )
        r_star = loader.get_r_star(V=2.0, eta=-1.0)
        phase_type = loader.get_phase_type(V=2.0, eta=-1.0)
    """

    def __init__(self, data: pd.DataFrame, theta_fixed: Optional[float] = None):
        """
        Initialize with Phase 1 data.

        Args:
            data: DataFrame with columns [V, eta, theta, r0, stability_type]
            theta_fixed: If specified, filter to this theta value
        """
        self.data = data
        self.theta_fixed = theta_fixed

        if theta_fixed is not None:
            theta_tol = 1e-6
            mask = np.abs(data['theta'] - theta_fixed) < theta_tol
            self.data = data[mask].copy()

            if len(self.data) == 0:
                raise ValueError(
                    f"No data for theta={theta_fixed}. "
                    f"Available: {sorted(data['theta'].unique())}"
                )

        stable_mask = (
            self.data['stability_type'].isin(['type-I', 'type-II']) &
            self.data['r0'].notna()
        )
        self.stable_data = self.data[stable_mask].copy()

        if len(self.stable_data) == 0:
            warnings.warn("No stable points in Phase 1 data.")
            self.r_star_interpolator = None
        else:
            self._build_interpolators()

    @classmethod
    def from_csv(cls, csv_path: Union[str, Path], theta_fixed: Optional[float] = None,
                 topology: Optional[str] = None):
        """Load Phase 1 results from CSV file.

        Args:
            csv_path: Path to CSV file
            theta_fixed: If specified, filter to this theta value
            topology: If specified, filter to this topology (e.g. 'S3', 'T3', 'Nil3').
                      Useful when loading a combined multi-topology CSV.
        """
        csv_path = Path(csv_path)
        if not csv_path.exists():
            raise FileNotFoundError(f"Phase 1 results not found: {csv_path}")

        data = pd.read_csv(csv_path)

        required_cols = ['V', 'eta', 'theta', 'r0', 'stability_type']
        missing = set(required_cols) - set(data.columns)
        if missing:
            raise ValueError(f"CSV missing columns: {missing}")

        if topology is not None and 'topology' in data.columns:
            data = data[data['topology'] == topology].copy()
            if len(data) == 0:
                raise ValueError(
                    f"No data for topology='{topology}' in {csv_path}. "
                    f"Available: {sorted(pd.read_csv(csv_path)['topology'].unique())}"
                )

        return cls(data, theta_fixed)

    def _build_interpolators(self):
        """Build 2D interpolators for r* and phase type."""
        eta_unique = np.sort(self.stable_data['eta'].unique())
        V_unique = np.sort(self.stable_data['V'].unique())

        expected_size = len(eta_unique) * len(V_unique)
        actual_size = len(self.stable_data)

        if actual_size == expected_size:
            self._build_regular_interpolators(eta_unique, V_unique)
        else:
            self._build_irregular_interpolators()

    def _build_regular_interpolators(self, eta_grid: np.ndarray, V_grid: np.ndarray):
        """Build interpolators for regular grid."""
        n_eta, n_V = len(eta_grid), len(V_grid)
        r_star_grid = np.full((n_eta, n_V), np.nan)
        phase_grid = np.full((n_eta, n_V), 0)

        self.phase_map = {'type-I': 1, 'type-II': 2, 'type-III': 3}
        self.phase_reverse_map = {1: 'I', 2: 'II', 3: 'III', 0: 'unknown'}

        for _, row in self.stable_data.iterrows():
            i = np.searchsorted(eta_grid, row['eta'])
            j = np.searchsorted(V_grid, row['V'])
            if i < n_eta and j < n_V:
                r_star_grid[i, j] = row['r0']
                phase_type = row['stability_type']
                phase_grid[i, j] = self.phase_map.get(phase_type, 0)

        self.r_star_interpolator = RegularGridInterpolator(
            (eta_grid, V_grid), r_star_grid,
            method='linear', bounds_error=False, fill_value=None
        )
        self.phase_interpolator = RegularGridInterpolator(
            (eta_grid, V_grid), phase_grid,
            method='nearest', bounds_error=False, fill_value=0
        )

        self.eta_bounds = (eta_grid[0], eta_grid[-1])
        self.V_bounds = (V_grid[0], V_grid[-1])
        self.grid_type = 'regular'

    def _build_irregular_interpolators(self):
        """Build interpolators for irregular grid."""
        points = self.stable_data[['eta', 'V']].values
        r_values = self.stable_data['r0'].values

        phase_map = {'type-I': 1, 'type-II': 2, 'type-III': 3}
        phase_values = self.stable_data['stability_type'].map(phase_map).values

        self.r_star_interpolator = LinearNDInterpolator(points, r_values, fill_value=np.nan)
        self.phase_interpolator = LinearNDInterpolator(points, phase_values, fill_value=0)

        self.eta_bounds = (points[:, 0].min(), points[:, 0].max())
        self.V_bounds = (points[:, 1].min(), points[:, 1].max())
        self.grid_type = 'irregular'
        self.phase_reverse_map = {1: 'I', 2: 'II', 3: 'III', 0: 'unknown'}

    def get_r_star(self, V: float, eta: float, warn_extrapolation: bool = True) -> Optional[float]:
        """Get stable radius r* at given (V, eta) point."""
        if self.r_star_interpolator is None:
            return None

        in_bounds = (
            self.eta_bounds[0] <= eta <= self.eta_bounds[1] and
            self.V_bounds[0] <= V <= self.V_bounds[1]
        )

        if not in_bounds and warn_extrapolation:
            warnings.warn(f"Query (V={V}, eta={eta}) outside bounds")

        if self.grid_type == 'regular':
            r_star = self.r_star_interpolator((eta, V))
        else:
            r_star = self.r_star_interpolator(eta, V)

        if isinstance(r_star, np.ndarray):
            r_star = r_star.item()

        if np.isnan(r_star):
            return None

        return float(r_star)

    def get_phase_type(self, V: float, eta: float) -> str:
        """Get phase type (I/II/III) at given (V, eta) point."""
        if self.phase_interpolator is None:
            return "unknown"

        if self.grid_type == 'regular':
            phase_int_float = self.phase_interpolator((eta, V))
        else:
            phase_int_float = self.phase_interpolator(eta, V)

        if isinstance(phase_int_float, np.ndarray):
            phase_int_float = phase_int_float.item()

        phase_int = int(round(float(phase_int_float)))
        return self.phase_reverse_map.get(phase_int, 'unknown')

    def is_in_stable_bounds(self, V: float, eta: float) -> bool:
        """Check if (V, eta) is within stable data bounds."""
        if not hasattr(self, 'eta_bounds') or not hasattr(self, 'V_bounds'):
            return False
        return (
            self.eta_bounds[0] <= eta <= self.eta_bounds[1] and
            self.V_bounds[0] <= V <= self.V_bounds[1]
        )

    def summary(self) -> Dict:
        """Get summary statistics of loaded Phase 1 data."""
        total = len(self.data)
        stable = len(self.stable_data)

        type_counts = self.data['stability_type'].value_counts().to_dict()

        return {
            'total_points': total,
            'stable_points': stable,
            'stable_fraction': stable / total if total > 0 else 0,
            'type_counts': type_counts,
            'eta_range': (self.data['eta'].min(), self.data['eta'].max()),
            'V_range': (self.data['V'].min(), self.data['V'].max()),
            'theta_values': sorted(self.data['theta'].unique()),
            'theta_fixed': self.theta_fixed,
            'grid_type': getattr(self, 'grid_type', 'none')
        }
