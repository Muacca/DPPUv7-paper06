# Curvature Layer

-> [English](README.md)

Curvature tensor computation、Pontryagin diagnostics、Hodge operations、Weyl tensors、および spatial-curvature helpers のためのモジュール群です。

## 概要

curvature layer は EC+NY engine と downstream callers が利用する symbolic tensor routines を提供します。この layer には次の要素が含まれます。

- `spatial_lie.py`: scale-stripped 3D scalar curvature と `chi` を導出します。
- `pontryagin_lz.py`: `P_form` check に使う Lorentzian form-Hodge Pontryagin density を計算します。

これらのモジュールは、form-Hodge Pontryagin density と script-level diagnostic contraction を分離して扱うための基盤です。

## モジュール

### `riemann.py`

Non-coordinate frame における Riemann tensor computation。strict antisymmetry checks と first-index lowering を含みます。

**主なクラス / 関数:**

- `RiemannAntisymmetryError`: strict verification が失敗した場合の exception。
- `verify_antisymmetry_strict(R, dim, ...)`: Riemann antisymmetry identities を検証します。
- `compute_riemann_tensor(gamma, structure_constants, dim, ...)`: `R^a_{bcd}` を計算します。
- `lower_first_index(Riemann, metric, dim)`: first index を下げて `R_abcd` を得ます。

### `ricci.py`

Ricci tensor、Ricci scalar、decomposition helpers。

**主な関数:**

- `compute_ricci_tensor(Riemann, dim)`: Ricci contraction。
- `compute_ricci_scalar(Riemann, dim, logger=None, metric_inv=None)`: Riemann data から Ricci scalar を計算します。
- `compute_ricci_scalar_from_tensor(Ricci, metric_inv, dim)`: Ricci tensor から scalar を計算します。
- `decompose_ricci_tensor(Ricci, metric, dim)`: trace/traceless decomposition。

### `hodge.py`

Hodge dual と 2-form block helpers。

**主な関数:**

- `hodge_dual_2form(F_down, metric, metric_inv=None, signature="lorentzian")`: 2-form の Hodge dual。
- `basis_2form(a, b)`: basis 2-form。
- `compute_hodge_dual(R_tensor)`: legacy Riemann Hodge dual helper。
- `classify_lz_2form_block(c, d)`: LZ-native block classification。
- `hodge_swaps_blocks()`: block-swap sanity check。

### `pontryagin.py`

Legacy numerical Pontryagin helper。

**主な関数:**

- `compute_P_from_riemann(R_arr)`: numerical array に対して `P=<R,*R>` を計算します。
- `get_riemann_numeric(...)`: lambdified engine result から numerical Riemann data を抽出します。

### `pontryagin_lz.py`

Lorentzian/LZ-native Pontryagin density helper。

**主な関数:**

- `pontryagin_density_lorentzian(R_abcd, metric, simplify_result=True)`: form-Hodge `P_form=<R,*R>` density を計算します。
- `pontryagin_density_lorentzian_by_epsilon(...)`: epsilon-contraction による cross-check route。
- `classify_pontryagin_expression(expr, zero_tol=None)`: exact-zero/nonzero symbolic result を分類します。
- `classify_pontryagin_expression_deep(expr, zero_tol=None)`: より深い simplification と classification route。

このモジュールは form-Hodge `P_form=<R,*R>` computation と exact-zero classification のためのものです。個別 script が実装する custom internal-pair diagnostic contraction とは別の量です。

### `spatial_lie.py`

Scale-stripped spatial Lie-frame curvature helper。

**主な関数:**

- `spatial_structure_constants_3d(params=None)`: five-parameter LZ family の 3D structure constants を構成します。
- `scalar_curvature_q2_from_structure_constants(structure_constants)`: `R^(3) q^2` を計算します。
- `chi_from_lz_parameters(params)`: `chi = R^(3) q^2 / 6` を計算します。

### `sd_extension.py`

Engine 用の self-duality（SD/ASD）mixin。

**主なクラス:**

- `SDExtensionMixin`: topology engine に SD functionality を付与します。

### `sd_diagnostics.py`

Self-duality と Pontryagin diagnostic evaluation。

**主なクラス / 関数:**

- `CurvatureSDDiagnostics`: SD/ASD residuals と関連 curvature diagnostics を評価します。

### `weyl.py`

Weyl tensor と Weyl scalar の計算。

**主な関数:**

- `compute_weyl_tensor(R_abcd, Ricci, R_scalar, metric, dim)`: `C_abcd` を計算します。
- `compute_weyl_scalar(C_abcd, metric)`: `C^2` を計算します。

## 注記

- `P_form` は `pontryagin_lz.py` で計算される form-Hodge density です。
- Internal-pair diagnostic contraction は script-level quantity であり、true form-Hodge Pontryagin density と混同しないでください。
- `chi` は `spatial_lie.py` により LC/spatial-curvature sector から導出されます。

## 依存関係

- [engine](../engine/README_ja.md): LC/EC connection data と frame pipeline。
- [topology](../topology/README_ja.md): LZ-native structure constants。
- [utils](../utils/README_ja.md): epsilon tensors と symbolic utilities。
