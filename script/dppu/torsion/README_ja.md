# Torsion Layer

-> [English](README.md)

Torsion modes、torsion tensor ansatz construction、torsion scalar、Nieh-Yan variants を定義するモジュール群です。

## 概要

torsion layer は EC+NY engine と reduced-sector helpers が使う `AX`, `VT`, `MX` torsion branches を提供します。runtime constructor は `construct_torsion_tensor(...)` です。この関数は legacy Euclidean convention と Lorentzian `lz_native` convention の両方をサポートします。

## モジュール

### `mode.py`

Torsion mode definitions。

| Mode | 説明 | Parameters |
|------|------|------------|
| `AX` | axial torsion only | `eta` |
| `VT` | vector-trace torsion only | `V` |
| `MX` | mixed torsion | `eta`, `V` |

```python
from dppu.torsion.mode import Mode

mode = Mode.MX
```

### `nieh_yan.py`

Nieh-Yan variant definitions と density computation。

| Variant | Definition | Meaning |
|---------|------------|---------|
| `TT` | `T^a wedge T_a` | torsion-torsion term |
| `REE` | `e^a wedge e^b wedge R_ab` | curvature-derived term |
| `FULL` | `TT - REE` | full Nieh-Yan density |

**主な API:**

- `NyVariant.TT`
- `NyVariant.REE`
- `NyVariant.FULL`
- `NyVariant.includes_tt`
- `NyVariant.includes_ree`

density 自体は engine pipeline の `E4.10` step で、利用可能な torsion tensor と curvature tensor から組み立てます。

### `ansatz.py`

Torsion tensor ansatz construction。

**主な関数:**

- `construct_torsion_tensor(mode, r, eta, V, metric, dim=4, logger=None, frame_convention="legacy_euclidean", signature="euclidean")`: `AX`, `VT`, `MX` に対応する `T^a_{bc}` を構成します。

Lorentzian LZ-native construction では次の形で使います。

```python
from sympy import Matrix, symbols
from dppu.torsion.ansatz import construct_torsion_tensor
from dppu.torsion.mode import Mode

q = symbols("q", positive=True)
eta, V = symbols("eta V", real=True)
metric = Matrix.diag(-1, 1, 1, 1)

T = construct_torsion_tensor(
    Mode.MX,
    q,
    eta,
    V,
    metric,
    dim=4,
    frame_convention="lz_native",
    signature="lorentzian",
)
```

### `scalar.py`

Torsion scalar computation。

```text
T = T_abc T^abc
```

## 注記

- `AX` と `VT` は single-channel active torsion branches です。
- `MX` は `eta` と `V` の両方を持ちます。
- Lorentzian LZ-native callers では `signature="lorentzian"` と `frame_convention="lz_native"` を使います。

## 依存関係

- [utils](../utils/README_ja.md): epsilon tensors と symbolic helpers。
- [engine](../engine/README_ja.md): contortion と EC connection construction。
- [curvature](../curvature/README_ja.md): `REE` Nieh-Yan term に使う curvature data。
