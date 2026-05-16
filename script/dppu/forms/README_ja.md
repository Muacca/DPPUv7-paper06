# Forms Layer

-> [English](README.md)

Exterior differential form algebra と Nieh-Yan form-engine helpers を提供します。

## 概要

forms layer は、軽量な symbolic exterior algebra representation と、自己完結した Nieh-Yan form engine を提供します。component-tensor engine pipeline とは分離されており、この package では ordered basis-index tuples の辞書として form を扱います。

## Form Representation

form は次の形で表現します。

```python
Form = dict[tuple[int, ...], object]
```

例:

- `{(): c}` は scalar 0-form。
- `{(0,): 1}` は `e^0`。
- `{(0, 1): c}` は `c e^0 wedge e^1`。

index tuple は昇順で保持され、wedge product の符号は coefficient に吸収されます。

## モジュール

### `exterior_algebra.py`

Pure exterior algebra operations。

**主な object / 関数:**

- `Form`: dictionary representation の type alias。
- `basis(i)`: elementary 1-form `e^i` を返します。
- `scalar_form(value)`: scalar を 0-form として包みます。
- `clean_form(form)`: simplification 後に zero-coefficient entries を除去します。
- `add_forms(*forms)`: forms を加算します。
- `scale_form(value, form)`: form に scalar を掛けます。
- `wedge(left, right)`: exterior product を計算します。
- `form_to_str(form)`: human-readable representation を返します。
- `coefficient(form, key=(0, 1, 2, 3))`: coefficient を抽出します。

### `nieh_yan.py`

Nieh-Yan boundary forms と exact-form checks のための exterior-form engine。

**主なクラス:**

- `NiehYanFormEngine`: symbol set を保持し、torsion forms、Nieh-Yan 3-form `Q`、`dQ`、`T^a wedge T_a`、`e^a wedge e^b wedge R_ab`、関連 audit data を計算します。

**主なメソッド:**

- `structure_constants()`: engine の coframe convention に対応する frame structure constants を構成します。
- `exterior_derivative(form, C=None)`: exterior derivative を計算します。
- `torsion_tensor_low(mode)`: `AX`, `VT`, `MX` に対応する lowered torsion tensor を構成します。
- `torsion_forms_low(mode)` / `torsion_forms_up(mode)`: torsion 2-forms を構成します。
- `q_nieh_yan(mode)`: `Q = e^a wedge T_a` を計算します。
- `integrated_q_boundary(mode)`: integrated boundary functional を計算します。
- `dq_reduced_density(mode)`: `dQ` から reduced density を計算します。
- `tt_form(mode)`: `T^a wedge T_a` を計算します。
- `curvature_forms_low(mode)`: EC curvature 2-forms を計算します。
- `e_e_r_form(mode)`: `e^a wedge e^b wedge R_ab` を計算します。
- `ny_form_audit(mode)`: `dQ` と `TT - e e R` を比較する dictionary を返します。

## 使用例

```python
from sympy import symbols
from types import SimpleNamespace

from dppu.forms import NiehYanFormEngine, basis, form_to_str, wedge

e01 = wedge(basis(0), basis(1))
print(form_to_str(e01))

r, rd, rdd, N, Nd, eta, etad, V, Vd, V3, kappa, theta_NY = symbols(
    "r rd rdd N Nd eta etad V Vd V3 kappa theta_NY"
)
sy = SimpleNamespace(
    r=r, rd=rd, rdd=rdd, N=N, Nd=Nd,
    eta=eta, etad=etad, V=V, Vd=Vd,
    V3=V3, kappa=kappa, theta_NY=theta_NY,
)

engine = NiehYanFormEngine(sy)
audit = engine.ny_form_audit("MX")
```

## `dppu.torsion.nieh_yan` との関係

Nieh-Yan 関連モジュールが 2 つあるのは意図的です。

- `dppu.torsion.nieh_yan` は engine pipeline において `TT`, `REE`, `FULL` の component-density variant を選択するための `NyVariant` enum を定義します。
- `dppu.forms.nieh_yan` は `Q`, `dQ`, exact-form diagnostics のための exterior-form engine である `NiehYanFormEngine` を定義します。

前者は variant-selection API、後者は form-calculation engine です。

## 依存関係

- SymPy: symbolic expressions and simplification。
- [utils](../utils/README_ja.md): epsilon symbols。
- [torsion](../torsion/README_ja.md): torsion mode naming convention。
