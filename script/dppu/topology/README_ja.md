# Topology Layer

-> [English](README.md)

Topology-specific computation engines と LZ-native topology invariants を提供するモジュール群です。

## 概要

topology layer は DPPUv7 EC+NY pipeline のために `M^3 x S^1` 型の幾何を実装します。runtime code では、明示的な topology、signature、frame convention、deformation、fiber、torsion fields を持つ unified configurable interface（`UnifiedEngine` + `DOFConfig`）を使います。

unified interface で利用できる canonical topology labels は次の 4 つです。

```text
S3 / T3 / Nil3 / Sol3
```

`Nil3` と `Sol3` はこの package では local Lie-frame reductions をサポートします。global quotient、spin-structure、spectral data は別の解析層です。

## クラス階層

```text
BaseFrameEngine            (dppu/engine/pipeline.py)
`-- TopologyEngine         (base_topology.py)
    |-- S3Geometry         (s3.py)
    |-- T3Geometry         (t3.py)
    |-- Nil3Geometry       (nil3.py)
    `-- Sol3Geometry       (sol3.py)
```

`UnifiedEngine` は `DOFConfig` に基づいて concrete implementation を選択します。

---

## モジュール

### `unified.py` - UnifiedEngine and DOFConfig

主要な runtime entry point です。

**主なクラス / 関数:**

- `TopologyType`: `S3`, `T3`, `NIL3`, `SOL3` を持つ enum。
- `FiberMode`: `NONE`, `TWIST`, `MIXING`, `BOTH` を持つ enum。
- `CurvatureSource`: `LC`, `EC` を持つ enum。
- `DOFConfig`: explicit configuration dataclass。
- `make_config(...)`: topology/configuration fields から `DOFConfig` を構成します。
- `make_engine(...)`: concrete topology engine を構成します。
- `normalize_topology_str(label)`: label を `S3`, `T3`, `Nil3`, `Sol3` に正規化します。

**Lorentzian LZ-native usage:**

Lorentzian LZ-native computation では `signature="lorentzian"` と `frame_convention="lz_native"` を使います。

```python
from dppu.topology import FiberMode, make_engine
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

engine = make_engine(
    "Sol3",
    torsion_mode=Mode.MX,
    ny_variant=NyVariant.FULL,
    signature="lorentzian",
    frame_convention="lz_native",
    fiber_mode=FiberMode.NONE,
    enable_squash=False,
)
engine.run()
```

### `base_topology.py` - TopologyEngine

Engine pipeline の topology-specific hooks を実装する abstract base class です。

**Subclass hooks:**

- `_build_radial_and_deformation_params(params)`: scale/deformation symbols を定義します。
- `_build_structure_constants(params, C)`: frame structure constants を埋めます。
- `_compute_volume(params)`: topology volume factor を返します。

### `s3.py` - S3Geometry

SU(2) Lie group に基づく `S3 x S1` engine。

- LZ reduced-sector specialization: `(a,b,c,u,v) = (4,4,4,0,0)`。
- LZ scalar `chi`: `4`。
- Closed isotropic reference topology。

### `t3.py` - T3Geometry

abelian flat three-torus に対応する `T3 x S1` engine。

- LZ reduced-sector specialization: `(a,b,c,u,v) = (0,0,0,0,0)`。
- LZ scalar `chi`: `0`。
- Flat reference topology。

### `nil3.py` - Nil3Geometry

Heisenberg nilmanifold に対応する `Nil3 x S1` engine。

- LZ reduced-sector specialization: `(a,b,c,u,v) = (0,0,-1,0,0)`。
- LZ scalar `chi`: `-1/12`。
- local isotropic-scale reduction をサポートします。
- `Nil3` は bi-invariant ではないため、general Koszul formula を使います。

### `sol3.py` - Sol3Geometry

solvable Thurston geometry に対応する `Sol3 x S1` engine。

- LZ reduced-sector specialization: `(a,b,c,u,v) = (0,0,0,1,-1)`。
- LZ scalar `chi`: `-1/3`。
- local isotropic-scale reduction をサポートします。
- ここで用いる diagonal frame scaling の下では spatial structure constants が rigid です。

### `lz_invariants.py`

LZ-native topology invariant extraction。

**主な object / 関数:**

- `A, B, C, U, W`: five-parameter LZ-native symbols。
- `Q`: positive isotropic scale symbol。
- `CANONICAL_TOPOLOGIES`: `("S3", "T3", "Nil3", "Sol3")`。
- `engine_lz_structure_constants(topology)`: LZ-native dense structure constants と engine scale symbol を返します。
- `extract_lz_five_parameter_specialization(topology)`: `(a,b,c,u,v)` を抽出します。
- `topology_lz_parameter_table()`: canonical specialization 全体を返します。

このモジュールは chi/scalar-curvature derivation と reduced-sector helpers から利用されます。

### `lz_adapter.py`

過去の convention transition のための legacy-to-LZ index conversion helpers と verification utilities。

### `bianchi_ix.py`

比較・guardrail 計算のための補助的な Bianchi IX / anisotropy helper。main `UnifiedEngine` topology dispatch には含まれません。

---

## Topology Comparison

| Topology | LZ specialization `(a,b,c,u,v)` | `chi` | note |
|----------|---------------------------------|------:|----------------|
| `S3` | `(4,4,4,0,0)` | `4` | closed isotropic reference |
| `T3` | `(0,0,0,0,0)` | `0` | flat reference |
| `Nil3` | `(0,0,-1,0,0)` | `-1/12` | local isotropic-scale reduction |
| `Sol3` | `(0,0,0,1,-1)` | `-1/3` | local isotropic-scale reduction |

## 依存関係

- [engine](../engine/README_ja.md): `BaseFrameEngine`, metric, connection, pipeline infrastructure。
- [curvature](../curvature/README_ja.md): curvature computation と spatial curvature helpers。
- [torsion](../torsion/README_ja.md): torsion mode construction。
- [action](../action/README_ja.md): action と effective-potential helpers。
