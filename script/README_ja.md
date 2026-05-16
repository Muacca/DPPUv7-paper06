# DPPUv7 paper06 スクリプトディレクトリ

-> [English](README.md)

**論文**: "Reduced-Sector chi-Universality in Lorentzian EC+NY Minisuperspace"（paper06）

Lorentzian Einstein-Cartan + Nieh-Yan reduced sector を、4 つの topology representative `S3`, `T3`, `Nil3`, `Sol3` 上で解析するための Python パッケージ群と実行スクリプトです。paper06 本文および Appendix B の再現性確認に対応し、admissibility classification、P-channel diagnostics、chi-bridge identity、reduced orbit atlas を支えます。

---

## ディレクトリ構成

```text
script/
|-- docs/                      # 技術ドキュメントと規約
|-- dppu/                      # メイン Python パッケージ（DPPUv7）
|   |-- action/                # EC action, Lagrangian, potential, paper06 reduced-sector helper
|   |-- analysis/              # Wick 回転分類と signature 解析
|   |-- curvature/             # Riemann, Ricci, Pontryagin, Weyl, SD, spatial Lie-frame curvature
|   |-- engine/                # metric, LC/EC connection, contortion, spin-2, pipeline
|   |-- forms/                 # exterior algebra と Nieh-Yan form
|   |-- kk/                    # Kaluza-Klein photon effective theory
|   |-- scanning/              # parameter scan, stability scan, SD diagnostics loader
|   |-- topology/              # topology engine と LZ-native topology invariants
|   |-- torsion/               # torsion mode, Ansatz, Nieh-Yan density
|   `-- utils/                 # 共通 utility と tee logging
`-- scripts/                   # 実行スクリプト
    |-- paper06/               # paper06 検証スクリプト
    `-- visualize/             # 図表生成スクリプト
```

### `docs/` - ドキュメント

技術ドキュメントと規約:

- [DPPUv7 Engine CONVENTIONS](docs/CONVENTIONS_ja.md) - engine core の規約と仕様
- [DPPUv7 SymPy guideline](docs/SymPy_guideline_ja.md) - SymPy 利用方針と実装ガイドライン
- [DPPUv7 LOGGING](docs/LOGGING_ja.md) - logging 規約と log format

---

## パッケージ概要（`dppu/`）

| モジュール | 役割 | 主なクラス / 関数 |
|------------|------|-------------------|
| [`action`](dppu/action/README_ja.md) | EC action, Lagrangian, effective potential, paper06 reduced-sector extraction | `reduced_sector`, `derive_static_branch_function`, `solve_auxiliary_shell` |
| [`analysis`](dppu/analysis/README_ja.md) | Wick/signature classification と diagnostic-comparison helpers | `wick_classification`, `build_phase1e_classification_table` |
| [`curvature`](dppu/curvature/README_ja.md) | curvature tensor, Pontryagin diagnostics, spatial Lie-frame curvature | `RiemannTensor`, `pontryagin`, `spatial_lie`, `sd_diagnostics` |
| [`engine`](dppu/engine/README_ja.md) | metric, LC/EC connection, contortion, spin-2, pipeline | `metric`, `levi_civita`, `ec_connection`, `contortion`, `pipeline` |
| [`forms`](dppu/forms/README_ja.md) | exterior differential form algebra と Nieh-Yan form engine | `Form`, `wedge`, `NiehYanFormEngine` |
| [`kk`](dppu/kk/README_ja.md) | KK photon effective theory | `extract_maxwell`, `extract_mass`, `extract_cs` |
| [`scanning`](dppu/scanning/README_ja.md) | parameter scan, stability diagnostics, result loader | `parameter_scan`, `stability`, `scan_results_loader` |
| [`topology`](dppu/topology/README_ja.md) | topology-specific engine と LZ five-parameter specialization | `UnifiedEngine`, `make_engine`, `lz_invariants` |
| [`torsion`](dppu/torsion/README_ja.md) | torsion mode と Nieh-Yan variant selection | `Mode`, `NyVariant`, `construct_torsion_tensor` |
| [`utils`](dppu/utils/README_ja.md) | 共通 utility と tee logging | `setup_log`, `teardown_log`, `epsilon_symbol`, `prove_zero` |

---

## 実行スクリプト概要（`scripts/`）

### `paper06/` - paper06 検証スクリプト

| スクリプト | 説明 |
|------------|------|
| `admissibility_classification.py` | `S3/T3/Nil3/Sol3 x EH/AX/VT/MX` の 4-topology branch classification を導出します。EH/AX/VT が `L_ADMISSIBLE`、MX が `L_CONDITIONALLY_ADMISSIBLE` であり、全 branch が `overall_status=PASS` になることを確認します。 |
| `pform_cancellation.py` | AX/VT/MX について、4 topology 上の Lorentzian form-Hodge Pontryagin density `P_form=<R,*R>` を計算し、block orthogonality による exact cancellation を確認します。 |
| `chi_bridge_symbolic.py` | MX internal-pair diagnostic の five-parameter LZ-native family identity を導出します。raw `P_int^MX` から `C` を抽出し、LC spatial curvature から独立に `chi` を導出して `C + 9*chi = 0` を確認します。 |
| `orbit_atlas.py` | auxiliary-shell reduced orbit atlas を導出します。`qdot^2 + chi = 0` を確認し、chi の符号による orbit sheet 分類、および reduced vacuum atlas 内に `BOUNCE_LIKE` / `RECOLLAPSE_LIKE` が現れないことを確認します。 |

#### `paper06/` 実行仕様

各 paper06 検証スクリプトは、`script/` から実行した場合に `../data/<script_name>_YYYYMMDD_HHMMSS.log`、すなわち paper06 直下の `data/` に log を生成します。現時点では command-line option は不要です。

| スクリプト | `script/` からの実行コマンド | 既定出力 | 期待される verdict |
|------------|------------------------------|----------|--------------------|
| `admissibility_classification.py` | `python scripts/paper06/admissibility_classification.py` | `../data/admissibility_classification_<timestamp>.log` | `overall_verdict=PASS` |
| `pform_cancellation.py` | `python scripts/paper06/pform_cancellation.py` | `../data/pform_cancellation_<timestamp>.log` | `overall_verdict=PASS` |
| `chi_bridge_symbolic.py` | `python scripts/paper06/chi_bridge_symbolic.py` | `../data/chi_bridge_symbolic_<timestamp>.log` | `overall_verdict=PASS` |
| `orbit_atlas.py` | `python scripts/paper06/orbit_atlas.py` | `../data/orbit_atlas_<timestamp>.log` | `overall_verdict=PASS` |

---

### `visualize/` - 図表

図生成スクリプトは実行ログを作成しません。各スクリプトは `../LaTeX/figures/` 配下の固定パスに `dpi=300` で出力し、現時点では command-line option を取りません。

| ファイル | 説明 | 出力先 |
|----------|------|--------|
| `fig01_three_level_architecture.py` | chi-controlled topology universality の three-level architecture：five-parameter family → scalar chi → admissibility / `P_int^MX` / reduced orbit（Figure 1、Section 1.2） | `../LaTeX/figures/fig01_three_level_architecture.png` |
| `fig02_chi_bridge_line.py` | Diagnostic bridge：4 topology の specialization を載せた `C_topology = -9*chi` 直線（Figure 2、Section 5.1） | `../LaTeX/figures/fig02_chi_bridge_line.png` |
| `fig03_loci_diagram.py` | Five-parameter family `(a,b,c,u,v)` における distinguished loci：class-A 部分空間 `(a,c)` と class-B 部分空間 `(u,v)`（Figure 3、Section 5.3） | `../LaTeX/figures/fig03_loci_diagram.png` |
| `fig04_chi_sign_orbit_atlas.py` | chi-sign reduced orbit atlas：chi > 0 / chi = 0 / chi < 0 の `q`-`t` 定性的 3 パネル模式図（Figure 4、Section 7.2） | `../LaTeX/figures/fig04_chi_sign_orbit_atlas.png` |

---

## クイックスタート

```bash
# paper06 の script ディレクトリへ移動
cd 40_paper/20_DPPUv7-paper06/script

# 依存パッケージをインストール
pip install -r requirements.txt

# 検証スクリプトを 1 本実行
python scripts/paper06/chi_bridge_symbolic.py

# bash で paper06 検証スクリプトをすべて実行
for f in scripts/paper06/*.py; do python "$f"; done

# paper06 の図をすべて再生成
for f in scripts/visualize/*.py; do python "$f"; done
```

PowerShell:

```powershell
Set-Location 40_paper/20_DPPUv7-paper06/script
pip install -r requirements.txt
Get-ChildItem scripts/paper06/*.py | Sort-Object Name | ForEach-Object { python $_.FullName }
```

---

## paper06 再現性マップ

| 論文内の位置 | 対応スクリプト |
|--------------|----------------|
| Section 3, Section 6: Hamiltonian and local admissibility taxonomy | `scripts/paper06/admissibility_classification.py` |
| Section 4.1: `P_form` exact cancellation | `scripts/paper06/pform_cancellation.py` |
| Section 4.3, Section 5: `P_int^MX` chi-specialisation and `C_topology=-9 chi` | `scripts/paper06/chi_bridge_symbolic.py` |
| Section 7, Appendix C: reduced orbit atlas and chi-sign separation | `scripts/paper06/orbit_atlas.py` |
| Appendix B: tables, computational checks, and reproducibility notes | all `scripts/paper06/*.py` logs |

## スコープ注記

- `P_int` は internal-pair diagnostic であり、true Pontryagin density ではありません。
- `P_form=0`, `H_PASS`, `L_ADMISSIBLE` は orbit stability の主張ではありません。
- `C_topology=-9 chi` と `qdot^2+chi=0` は reduced-sector result であり、full-theory topology-universality theorem ではありません。
- `Nil3` と `Sol3` の entries は paper06 で用いる isotropic-scale reduction を指します。anisotropic extension と global quotient extension は将来課題です。
