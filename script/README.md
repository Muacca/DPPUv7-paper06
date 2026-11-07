# Computational supplement

[日本語](README_ja.md) · [Conventions](CONVENTIONS.md) · [Earlier conventions](legacy_conventions.md)

Use `paper06.py recompute --output-dir NEW` to derive the exact results and figure data. Use `paper06.py render --data-dir DATA --output-dir NEW` to draw figures from existing data. Run from the repository root with the `script/` prefix, as shown in the main README. Both commands require a new destination directory and fail rather than overwrite it. Install `requirements.txt` with Python 3.11 or later.

The public modules are:

| Module | Content |
| --- | --- |
| `geometry` | Sparse exterior forms; class-A connection; torsion, curvature, NY, CS and transgression |
| `derive` | Symbolic identities, bulk/endpoint reduction, Euler equations, constraints and four-geometry specialization |
| `shape`, `torus` | General shape variation, mixed cone and exact torus family |
| `invariants` | Quadratic/cubic contractions and Ricci inverse problem |
| `evolution` | Diagonal vacuum ODE, Taylor jets and 12 numerical points |
| `frames`, `canonical` | Global frames and endpoint canonical maps |
| `components`, `identities` | Koszul calculation and identity checks |
| `results` | Recompute and serialize exact data, dynamic regimes, field equations and figure curves |
| `render` | Read saved JSON/CSV and render six PNG/PDF figure pairs |

`geometry.class_a()` returns forms and the all-lowered tensor `R[a,b,c,d]`. `derive.derive_geometry()` returns this geometry and the derived scalar/endpoint expressions. `derive.derive_reduction(geometry_data)` derives the coordinate-time action. `derive.derive_classification(geometry_data)` specializes the fixed spatial metrics. The symbolic inputs and scopes are defined in the conventions guide. These functions do not classify general anisotropic solutions or physical spectra.

Saved-data rendering does not run the derivation. Symbolic recomputation constructs the connection before restricting profiles; it does not merely substitute numbers into the manuscript formulas. The supplement runs with the supplied files and the listed Python dependencies.

## Environment and numerical scope

Tested: Python 3.12.10; SymPy 1.14.0; NumPy 2.3.5; SciPy 1.16.3; Matplotlib 3.10.7. The pinned dependencies require Python 3.11 or later. Set `SYMPY_GROUND_TYPES=python` before starting Python to reproduce the tested symbolic backend (`$env:SYMPY_GROUND_TYPES='python'` in PowerShell; `export SYMPY_GROUND_TYPES=python` in a POSIX shell).

`recompute` constructs the class-A connection and the independent torsionless shape action, obtains contractions from fresh curvature and integrates the diagonal ODE with SciPy DOP853. The 12 points use q0=1, both H0 signs for Nil and Sol, times 0.01, 0.02, 0.04, rtol=2e-13, atol=2e-15, and max_step=0.002. These are point checks, not certified interval bounds. Symbolic expressions are exact; floating-point output may vary slightly with platforms.

`render` reads JSON/CSV and only evaluates the stored expression strings numerically; it runs no curvature derivation or ODE integration. The arguments of the main text and appendices do not depend on running the code.

The three folders are self-contained. Both commands refuse an existing output directory.
