# Scanning Layer

⇒ [日本語](README_ja.md)

Module group for parameter space scanning and scan-result integration.

## Overview

Provides (V, eta, theta) parameter space scanning, scan result loading/interpolation, and SD diagnostics integration.

This layer is currently wired for `S3/T3/Nil3` scans. `Sol3` support exists in the topology layer but is not integrated into these scan helpers.

## Modules

### parameter_scan.py

Grid scanning of parameter space.

**Key functions:**

- `run_scan(...)`: Scan parameter space

**Usage:**

```python
from dppu.scanning import run_scan

results = run_scan(
    V_points=50, eta_points=100, theta_points=20,
    V_min=0.0, V_max=5.0,
    eta_min=-5.0, eta_max=5.0,
    theta_min=0.0, theta_max=5.0,
    topologies=['S3', 'T3', 'Nil3'],
    ny_variants=['FULL', 'TT', 'REE'],
    output_dir='output/',
    n_workers=8
)
```

**Output CSV format:**

| Column | Description |
|--------|-------------|
| topology | S3, T3, Nil3 |
| ny_variant | FULL, TT, REE |
| V, eta, theta | Parameter values |
| r0 | Stable point position (None if type-III) |
| delta_V | Barrier height or well depth |
| stability_type | type-I, type-II, type-III |

### scan_results_loader.py

Loading and interpolation of scan results.

**ScanResultsLoader:**

```python
from dppu.scanning import ScanResultsLoader

# Load from CSV
loader = ScanResultsLoader.from_csv(
    'output/dppu_scan_S3_FULL.csv',
    theta_fixed=0.0
)

# Display summary
print(loader.summary())
# {'total_points': 5000, 'stable_points': 3200, ...}

# Get interpolated r*
r_star = loader.get_r_star(V=2.0, eta=-1.0)

# Get stability type
phase = loader.get_phase_type(V=2.0, eta=-1.0)
# 'I', 'II', or 'III'
```

**Grid types:**

- `regular`: Regular grid → Uses RegularGridInterpolator
- `irregular`: Irregular grid → Uses LinearNDInterpolator

### scan_sd_diagnostics.py

SD diagnostics integrated with scan results.

**SDScanDiagnostics:**

```python
from dppu.scanning import SDScanDiagnostics, ScanResultsLoader

loader = ScanResultsLoader.from_csv('scan.csv', theta_fixed=0.0)
sd_diag = SDScanDiagnostics(engine, scan_loader=loader)

# Evaluate SD at r*
result = sd_diag.evaluate_at_rstar(V=2.0, eta=-1.0, theta_NY=0.0)
print(f"r* = {result.r_star:.3f}")
print(f"SD residual = {result.sd_residual:.4f}")

# Scan the (η, V) plane
results = sd_diag.scan_with_rstar(
    eta_range=(-5, 3, 100),
    V_range=(0.5, 5, 50),
    theta_NY=0.0
)

print(f"SD curve points: {len(results['sd_curve'])}")
print(f"Type I ∩ SD: {len(results['type_I_sd_intersection'])}")
```

**SDScanResult:**

```python
@dataclass
class SDScanResult:
    eta: float
    V: float
    r_star: float
    phase_type: str
    sd_residual: float
    asd_residual: float
    curvature_norm: float
    is_nontrivial_sd: bool
    is_nontrivial_asd: bool
```

### stability.py

Stability classification of effective potential minima.

**Stability types:**

| Type | Condition | Physical Interpretation |
|------|-----------|------------------------|
| type-I | V''(r*) > 0, V(r*) > 0 | Metastable (with barrier) |
| type-II | V''(r*) > 0, V(r*) < 0 | True minimum (below vacuum) |
| type-III | No minimum or V''(r*) ≤ 0 | Unstable |

**Key functions:**

- `classify_stability(V, r_star)`: Classify stability
- `StabilityResult`: Data class for classification results

```python
from dppu.scanning import classify_stability

result = classify_stability(V_eff, r_star=1.5)
print(f"Type: {result.stability_type}")
print(f"r* = {result.r_star:.3f}")
```

## Workflow

### Scan -> SD Diagnostics Integration

1. **Scan**: parameter scan to identify stable points `r*(V, eta)`.
2. **Diagnostics**: execute SD diagnostics at `r*`.
3. **Analysis**: Investigate intersection of SD curve and stable region

Use the Python APIs (`run_scan`, `ScanResultsLoader`, and `SDScanDiagnostics`) from a project-specific driver script.

## Dependencies

- [topology](../topology/README.md): Computation engines
- [curvature](../curvature/README.md): SD diagnostics
- pandas (CSV loading)
- scipy (interpolation)

## Related Modules

- [action](../action/README.md): Lagrangian and potential computation
