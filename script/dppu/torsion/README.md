# Torsion Layer

-> [Japanese](README_ja.md)

Module group defining torsion modes, torsion tensor ansatz construction, torsion scalar, and Nieh-Yan variants.

## Overview

The torsion layer supplies the `AX`, `VT`, and `MX` torsion branches used by the EC+NY engine and reduced-sector helpers.  The runtime constructor is `construct_torsion_tensor(...)`; it supports both the legacy Euclidean convention and the Lorentzian `lz_native` convention.

## Modules

### `mode.py`

Torsion mode definitions.

| Mode | Description | Parameters |
|------|-------------|------------|
| `AX` | axial torsion only | `eta` |
| `VT` | vector-trace torsion only | `V` |
| `MX` | mixed torsion | `eta`, `V` |

```python
from dppu.torsion.mode import Mode

mode = Mode.MX
```

### `nieh_yan.py`

Nieh-Yan variant definitions and density computation.

| Variant | Definition | Meaning |
|---------|------------|---------|
| `TT` | `T^a wedge T_a` | torsion-torsion term |
| `REE` | `e^a wedge e^b wedge R_ab` | curvature-derived term |
| `FULL` | `TT - REE` | full Nieh-Yan density |

**Key API:**

- `NyVariant.TT`
- `NyVariant.REE`
- `NyVariant.FULL`
- `NyVariant.includes_tt`
- `NyVariant.includes_ree`

The density itself is assembled in the engine pipeline step `E4.10` from the available torsion and curvature tensors.

### `ansatz.py`

Torsion tensor ansatz construction.

**Key function:**

- `construct_torsion_tensor(mode, r, eta, V, metric, dim=4, logger=None, frame_convention="legacy_euclidean", signature="euclidean")`: construct `T^a_{bc}` for `AX`, `VT`, or `MX`.

For Lorentzian LZ-native construction, use:

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

Torsion scalar computation.

```text
T = T_abc T^abc
```

## Notes

- `AX` and `VT` are the single-channel active torsion branches.
- `MX` carries both `eta` and `V`.
- Lorentzian LZ-native callers should use `signature="lorentzian"` and `frame_convention="lz_native"`.

## Dependencies

- [utils](../utils/README.md): epsilon tensors and symbolic helpers.
- [engine](../engine/README.md): contortion and EC connection construction.
- [curvature](../curvature/README.md): curvature data for the `REE` Nieh-Yan term.
