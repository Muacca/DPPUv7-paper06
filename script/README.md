# DPPUv7 paper06 Scripts Directory

-> [Japanese](README_ja.md)

**Paper**: "Reduced-Sector chi-Universality in Lorentzian EC+NY Minisuperspace" (paper06)

Python packages and execution scripts for the paper06 Lorentzian Einstein-Cartan + Nieh-Yan reduced-sector analysis over the four topology representatives `S3`, `T3`, `Nil3`, and `Sol3`.  The scripts support the admissibility classification, P-channel diagnostics, chi-bridge identity, and reduced orbit atlas used in the paper draft and Appendix B reproducibility notes.

---

## Directory Structure

```text
script/
|-- docs/                      # Technical documentation and conventions
|-- dppu/                      # Main Python package (DPPUv7)
|   |-- action/                # EC action, Lagrangian, potentials, paper06 reduced-sector helpers
|   |-- analysis/              # Wick classification and signature analysis
|   |-- curvature/             # Riemann, Ricci, Pontryagin, Weyl, SD, spatial Lie-frame curvature
|   |-- engine/                # Metric, LC/EC connection, contortion, spin-2 sector, pipeline
|   |-- forms/                 # Exterior algebra and Nieh-Yan form
|   |-- kk/                    # Kaluza-Klein photon effective theory
|   |-- scanning/              # Parameter scan, stability scan, SD diagnostics loader
|   |-- topology/              # Topology engines and LZ-native topology invariants
|   |-- torsion/               # Torsion modes, Ansatz, Nieh-Yan density
|   `-- utils/                 # Common utilities and tee logging
`-- scripts/                   # Execution scripts
    |-- paper06/               # paper06 verification scripts
    `-- visualize/             # Figure generation scripts
```

### `docs/` - Documentation

Technical documentation and conventions:

- [DPPUv7 Engine CONVENTIONS](docs/CONVENTIONS.md) - engine core conventions and specifications
- [DPPUv7 SymPy guideline](docs/SymPy_guideline.md) - SymPy usage guidelines and best practices
- [DPPUv7 LOGGING](docs/LOGGING.md) - logging conventions and log format specifications

---

## Package Overview (`dppu/`)

| Module | Role | Key Classes / Functions |
|--------|------|-------------------------|
| [`action`](dppu/action/README.md) | EC action, Lagrangian, effective potential, paper06 reduced-sector extraction | `reduced_sector`, `derive_static_branch_function`, `solve_auxiliary_shell` |
| [`analysis`](dppu/analysis/README.md) | Wick/signature classification and diagnostic-comparison helpers | `wick_classification`, `build_phase1e_classification_table` |
| [`curvature`](dppu/curvature/README.md) | Curvature tensors, Pontryagin diagnostics, spatial Lie-frame curvature | `RiemannTensor`, `pontryagin`, `spatial_lie`, `sd_diagnostics` |
| [`engine`](dppu/engine/README.md) | Metric, LC/EC connection, contortion, spin-2, pipeline | `metric`, `levi_civita`, `ec_connection`, `contortion`, `pipeline` |
| [`forms`](dppu/forms/README.md) | Exterior differential form algebra and Nieh-Yan form engine | `Form`, `wedge`, `NiehYanFormEngine` |
| [`kk`](dppu/kk/README.md) | KK photon effective theory | `extract_maxwell`, `extract_mass`, `extract_cs` |
| [`scanning`](dppu/scanning/README.md) | Parameter scan, stability diagnostics, result loader | `parameter_scan`, `stability`, `scan_results_loader` |
| [`topology`](dppu/topology/README.md) | Topology-specific engines and LZ five-parameter specializations | `UnifiedEngine`, `make_engine`, `lz_invariants` |
| [`torsion`](dppu/torsion/README.md) | Torsion modes and Nieh-Yan variant selection | `Mode`, `NyVariant`, `construct_torsion_tensor` |
| [`utils`](dppu/utils/README.md) | Common utilities and tee logging | `setup_log`, `teardown_log`, `epsilon_symbol`, `prove_zero` |

---

## Script Overview (`scripts/`)

### `paper06/` - paper06 Verification Scripts

| Script | Description |
|--------|-------------|
| `admissibility_classification.py` | Derives the four-topology branch classification for `S3/T3/Nil3/Sol3 x EH/AX/VT/MX`.  It verifies that EH/AX/VT are `L_ADMISSIBLE`, MX is `L_CONDITIONALLY_ADMISSIBLE`, and every branch has `overall_status=PASS`. |
| `pform_cancellation.py` | Computes the Lorentzian form-Hodge Pontryagin density `P_form=<R,*R>` for AX/VT/MX across the four topologies and checks exact cancellation by block orthogonality. |
| `chi_bridge_symbolic.py` | Derives the five-parameter LZ-native family identity for the MX internal-pair diagnostic.  It extracts `C` from raw `P_int^MX`, derives `chi` independently from LC spatial curvature, and checks `C + 9*chi = 0`. |
| `orbit_atlas.py` | Derives the auxiliary-shell reduced orbit atlas.  It verifies `qdot^2 + chi = 0` and classifies the chi-sign orbit sheets, including the absence of `BOUNCE_LIKE` and `RECOLLAPSE_LIKE` entries in this reduced vacuum atlas. |

#### `paper06/` Execution Matrix

Each paper06 verification script writes to `../data/<script_name>_YYYYMMDD_HHMMSS.log` when run from `script/` (that is, the `data/` directory at the paper06 root).  The scripts do not currently require command-line options.

| Script | Command from `script/` | Default output | Expected verdict |
|--------|------------------------|----------------|------------------|
| `admissibility_classification.py` | `python scripts/paper06/admissibility_classification.py` | `../data/admissibility_classification_<timestamp>.log` | `overall_verdict=PASS` |
| `pform_cancellation.py` | `python scripts/paper06/pform_cancellation.py` | `../data/pform_cancellation_<timestamp>.log` | `overall_verdict=PASS` |
| `chi_bridge_symbolic.py` | `python scripts/paper06/chi_bridge_symbolic.py` | `../data/chi_bridge_symbolic_<timestamp>.log` | `overall_verdict=PASS` |
| `orbit_atlas.py` | `python scripts/paper06/orbit_atlas.py` | `../data/orbit_atlas_<timestamp>.log` | `overall_verdict=PASS` |

---

### `visualize/` - Figures

Figure scripts do not create execution logs.  Each script writes to a fixed path under `../LaTeX/figures/` at `dpi=300` and does not currently take command-line options.

| File | Description | Output |
|------|-------------|--------|
| `fig01_three_level_architecture.py` | Three-level architecture of chi-controlled topology universality: five-parameter family → scalar chi → admissibility / `P_int^MX` / reduced orbit (Figure 1, Section 1.2) | `../LaTeX/figures/fig01_three_level_architecture.png` |
| `fig02_chi_bridge_line.py` | Diagnostic bridge: the `C_topology = -9*chi` line with the four topology specialisations (Figure 2, Section 5.1) | `../LaTeX/figures/fig02_chi_bridge_line.png` |
| `fig03_loci_diagram.py` | Distinguished loci in the five-parameter family `(a,b,c,u,v)`: class-A subspace `(a,c)` and class-B subspace `(u,v)` (Figure 3, Section 5.3) | `../LaTeX/figures/fig03_loci_diagram.png` |
| `fig04_chi_sign_orbit_atlas.py` | chi-sign reduced orbit atlas, three-panel qualitative `q`-`t` schematic for chi > 0, chi = 0, and chi < 0 (Figure 4, Section 7.2) | `../LaTeX/figures/fig04_chi_sign_orbit_atlas.png` |

---

## Quick Start

```bash
# From the paper06 script directory
cd 40_paper/20_DPPUv7-paper06/script

# Install dependencies
pip install -r requirements.txt

# Run one verification script
python scripts/paper06/chi_bridge_symbolic.py

# Run all paper06 verification scripts in bash
for f in scripts/paper06/*.py; do python "$f"; done

# Regenerate every paper06 figure
for f in scripts/visualize/*.py; do python "$f"; done
```

PowerShell:

```powershell
Set-Location 40_paper/20_DPPUv7-paper06/script
pip install -r requirements.txt
Get-ChildItem scripts/paper06/*.py | Sort-Object Name | ForEach-Object { python $_.FullName }
```

---

## Paper06 Reproducibility Map

| Paper location | Supporting script |
|----------------|-------------------|
| Section 3, Section 6: Hamiltonian and local admissibility taxonomy | `scripts/paper06/admissibility_classification.py` |
| Section 4.1: `P_form` exact cancellation | `scripts/paper06/pform_cancellation.py` |
| Section 4.3, Section 5: `P_int^MX` chi-specialisation and `C_topology=-9 chi` | `scripts/paper06/chi_bridge_symbolic.py` |
| Section 7, Appendix C: reduced orbit atlas and chi-sign separation | `scripts/paper06/orbit_atlas.py` |
| Appendix B: tables, computational checks, and reproducibility notes | all `scripts/paper06/*.py` logs |

## Scope Notes

- `P_int` is an internal-pair diagnostic, not a true Pontryagin density.
- `P_form=0`, `H_PASS`, and `L_ADMISSIBLE` are not orbit-stability claims.
- `C_topology=-9 chi` and `qdot^2+chi=0` are reduced-sector results, not full-theory topology-universality theorems.
- `Nil3` and `Sol3` entries refer to the isotropic-scale reductions used in paper06; anisotropic and global-quotient extensions remain future work.
