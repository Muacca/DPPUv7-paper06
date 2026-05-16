# KK レイヤー

⇒ [English](README.md)

Kaluza-Klein 光子有効理論のモジュール群。

## 概要

5次元 Einstein-Cartan 重力から O(A²) 光子有効作用を2つの独立ルート（高速ショートカットと完全 Riemann 検証）で計算する。全3トポロジーに対して Maxwell・Chern-Simons・質量² 係数の抽出を提供する。

**有効作用の構造**（内部多様体で積分後）:

```
S_eff ∝ ∫ d⁴x [ Σ_{i<j} k_ij F²_{ij}  +  k_CS ε^{ijk} A_i ∂_j A_k  +  Σ_i m²_i A_i² ]
```

---

## 2つの計算ルート

```
ショートカットルート（高速）:
    field_strength → gamma_gamma → extractor

完全検証ルート（低速）:
    field_strength → full_riemann → extractor
    └── ショートカット結果と相互検証
```

---

## モジュール

### field_strength.py

トポロジー対応の修正場の強さ F̃ と O(A¹) 接続摂動 ω⁽¹⁾。

**主要関数:**

- `make_F_plain(j, k, A, dA)`: 通常の場の強さ F_{jk} = ∂_j A_k − ∂_k A_j
- `make_F_tilde(j, k, A, dA, corrections)`: トポロジー固有補正を含む修正場の強さ F̃
- `s3_corrections(A)`: S³ 補正 dict（SU(2) 構造定数）
- `nil3_corrections(A)`: Nil³ 補正 dict（Heisenberg 構造定数）
- `t3_corrections(A)`: T³ 補正 dict（空 — 平坦、補正なし）
- `sol3_corrections(A, scale=1)`: homogeneous Sol 構造定数に対応する Sol³ 補正 dict
- `sol3_inhomogeneous_corrections(A, R_inv)`: 任意の `1/R` profile を持つ Sol³ 補正 dict
- `make_omega1(F_tilde_fn, r0, L)`: Koszul 公式による O(A¹) 接続摂動の構築
- `omega1_to_array(omega1_fn, n)`: ω⁽¹⁾ 関数を3次元配列に変換

### gamma_gamma.py

O(A²) Ricci スカラーの Γ×Γ ショートカット — 高速ルート。

**主要関数:**

- `gamma_gamma_ricci(omega1_fn)`: R^{Γ²} = Σ_{a,b,c} ω⁽¹⁾_{abc} ω⁽¹⁾_{abc} を計算（O(A²) で厳密）

**注記:** O(A²) において完全 Riemann ルートと代数的に等価で、約 21 倍高速。

### full_riemann.py

完全 Riemann 検証ルート — 低速だがショートカットの相互検証に使用。

**主要関数:**

- `make_omega0(topology, r0, L)`: 指定トポロジーの背景（O(A⁰)）接続の構築
- `full_riemann_scalar_a2(omega0, omega1_fn)`: 完全 R[ω⁽⁰⁾+ω⁽¹⁾] の O(A²) 成分を計算

**符号規約:** エンジン全体で使用する `[E_b, E_c] = −C^a_{bc} E_a` 規約に合わせ `−R_a2` を返す。

### extractor.py

Ricci スカラーからの Maxwell・Chern-Simons・質量² 係数抽出。

**主要関数:**

- `extract_maxwell(R_scalar, dA)`: 各 F²_{ij} 方向の `{(i,j): coeff}` を抽出
- `extract_mass(R_scalar, A, dA)`: 各 A_i² 質量項の `{i: coeff}` を抽出
- `extract_cs(R_scalar, A, dA)`: ε^{ijk} A_i ∂_j A_k 項の k_CS を抽出（なければ `None`）
- `extract_all(R_scalar, A, dA)`: `{'maxwell': ..., 'mass': ..., 'cs': ...}` を返すラッパー

**戻り値の規約:** 異方的な空間（squashed S³ など）を正確に扱うため、全抽出関数は dict を返す。空の dict はその項が存在しないことを意味する。

### validator.py

相互ルート検証ユーティリティ。

**主要関数:**

- `validate_kk_routes(topology)`: 指定トポロジーで両ルートを実行し、一致を検証

---

## 使用例

### クイックスタート（ショートカットルート）

```python
from sympy import symbols, Symbol
from dppu.kk import (
    s3_corrections, make_F_tilde, make_omega1,
    gamma_gamma_ricci, extract_all,
)

r0, L = symbols('r0 L', positive=True)
A   = [Symbol('A0'), Symbol('A1'), Symbol('A2')]
dA  = [[Symbol(f'dA{j}{k}') for k in range(3)] for j in range(3)]

# S³ トポロジー
corrections = s3_corrections(A)
F_tilde_fn  = lambda j, k: make_F_tilde(j, k, A, dA, corrections)
omega1_fn   = make_omega1(F_tilde_fn, r0, L)
R_GG        = gamma_gamma_ricci(omega1_fn)
coeffs      = extract_all(R_GG, A, dA)

print(coeffs['maxwell'])   # {(0,1): ..., (0,2): ..., (1,2): ...}
print(coeffs['mass'])      # {0: ..., 1: ..., 2: ...}  （S³：3重縮退）
print(coeffs['cs'])        # k_CS の式
```

### 検証

```python
from dppu.kk.validator import validate_kk_routes

validate_kk_routes('t3')    # T³ — 高速
validate_kk_routes('nil3')  # Nil³
validate_kk_routes('s3')    # S³ — 低速
```

`field_strength.py` には Sol³ 補正 helper も含まれます。ただし route validator は現時点では `S3/T3/Nil3` に限定されています。

---

## トポロジー別結果まとめ

| トポロジー | Maxwell k | CS k | 質量² |
|-----------|-----------|------|------|
| T³×S¹ | −L²/(2r₀⁴)（均一）| 0 | 0（質量なし）|
| Nil³×S¹ | −L²/(2r₀⁴)（均一）| −L²/r₀⁴ × A₂∂₀A₁（01 方向のみ）| −L²/(2r₀⁴)（A₂ のみ）|
| S³×S¹ | −L²/(2r₀⁴)（均一）| −2L²/r₀⁴（全3方向）| −2L²/r₀⁴（3重縮退）|

Maxwell 係数は全3トポロジーで同一の普遍値 −L²/(2r₀⁴) となる。

---

## 依存関係

- [engine](../engine/README_ja.md): ω⁽¹⁾ 用 Koszul 公式
- [topology](../topology/README_ja.md): 背景接続 ω⁽⁰⁾
- [utils](../utils/README_ja.md): エプシロン記号（`epsilon_3d`）

## 関連利用

- Shortcut route と full-Riemann route の比較には `dppu.kk.validator.validate_kk_routes(...)` を使います。
- Symbolic Ricci scalar から Maxwell、Chern-Simons、mass coefficients を抽出する script では `dppu.kk.extractor.extract_all(...)` を使います。
