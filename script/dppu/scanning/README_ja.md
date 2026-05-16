# Scanning Layer

⇒ [English](README.md)

パラメータ空間の走査と scan-result integration を担当するモジュール群。

## 概要

(V, eta, theta) パラメータ空間のスキャン、scan result の読み込み・補間、SD 診断との統合を提供。

この layer は現在の実装では `S3/T3/Nil3` scans に接続されています。`Sol3` support は topology layer にはありますが、この scan helper 群には未統合です。

## モジュール

### parameter_scan.py

パラメータ空間のグリッドスキャン。

**主要関数:**

- `run_scan(...)`: パラメータ空間をスキャン

**使用例:**

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

**出力CSV形式:**

| カラム | 説明 |
|--------|------|
| topology | S3, T3, Nil3 |
| ny_variant | FULL, TT, REE |
| V, eta, theta | パラメータ値 |
| r0 | 安定点の位置 (None if type-III) |
| delta_V | バリア高さまたは井戸深さ |
| stability_type | type-I, type-II, type-III |

### scan_results_loader.py

スキャン結果の読み込みと補間。

**ScanResultsLoader:**

```python
from dppu.scanning import ScanResultsLoader

# CSVから読み込み
loader = ScanResultsLoader.from_csv(
    'output/dppu_scan_S3_FULL.csv',
    theta_fixed=0.0
)

# サマリー表示
print(loader.summary())
# {'total_points': 5000, 'stable_points': 3200, ...}

# r*の補間取得
r_star = loader.get_r_star(V=2.0, eta=-1.0)

# 安定性タイプ取得
phase = loader.get_phase_type(V=2.0, eta=-1.0)
# 'I', 'II', or 'III'
```

**グリッドタイプ:**

- `regular`: 規則グリッド → RegularGridInterpolator使用
- `irregular`: 不規則グリッド → LinearNDInterpolator使用

### scan_sd_diagnostics.py

SD診断とスキャン結果の統合。

**SDScanDiagnostics:**

```python
from dppu.scanning import SDScanDiagnostics, ScanResultsLoader

loader = ScanResultsLoader.from_csv('scan.csv', theta_fixed=0.0)
sd_diag = SDScanDiagnostics(engine, scan_loader=loader)

# r*でのSD評価
result = sd_diag.evaluate_at_rstar(V=2.0, eta=-1.0, theta_NY=0.0)
print(f"r* = {result.r_star:.3f}")
print(f"SD residual = {result.sd_residual:.4f}")

# (η, V)平面のスキャン
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

有効ポテンシャル極小点の安定性分類。

**安定性タイプ:**

| タイプ | 条件 | 物理的解釈 |
|--------|------|-----------|
| type-I | V''(r*) > 0, V(r*) > 0 | 準安定（バリアあり） |
| type-II | V''(r*) > 0, V(r*) < 0 | 真の極小（真空より低い） |
| type-III | 極小なし or V''(r*) ≤ 0 | 不安定 |

**主要関数:**

- `classify_stability(V, r_star)`: 安定性を分類
- `StabilityResult`: 分類結果のデータクラス

```python
from dppu.scanning import classify_stability

result = classify_stability(V_eff, r_star=1.5)
print(f"Type: {result.stability_type}")
print(f"r* = {result.r_star:.3f}")
```

## ワークフロー

### Scan -> SD Diagnostics Integration

1. **Scan**: パラメータスキャンで安定点 `r*(V, eta)` を特定
2. **Diagnostics**: `r*` での SD 診断を実行
3. **分析**: SD曲線と安定領域の交差を調査

Project-specific driver script から Python APIs（`run_scan`, `ScanResultsLoader`, `SDScanDiagnostics`）を呼び出して使います。

## 依存関係

- [topology](../topology/README_ja.md): 計算エンジン
- [curvature](../curvature/README_ja.md): SD診断
- pandas (CSV読み込み)
- scipy (補間)

## 関連モジュール

- [action](../action/README_ja.md): ラグランジアンとポテンシャル計算
