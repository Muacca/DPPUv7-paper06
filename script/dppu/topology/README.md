# Topology Layer

-> [Japanese](README_ja.md)

Module group providing topology-specific computation engines and LZ-native topology invariants.

## Overview

The topology layer implements the geometry of `M^3 x S^1` manifolds for the DPPUv7 EC+NY pipeline.  Runtime code uses the unified configurable interface (`UnifiedEngine` + `DOFConfig`) with explicit topology, signature, frame convention, deformation, fiber, and torsion fields.

The canonical topology labels available through the unified interface are:

```text
S3 / T3 / Nil3 / Sol3
```

`Nil3` and `Sol3` support local Lie-frame reductions in this package.  Global quotient, spin-structure, and spectral data are separate layers of analysis.

## Class Hierarchy

```text
BaseFrameEngine            (dppu/engine/pipeline.py)
`-- TopologyEngine         (base_topology.py)
    |-- S3Geometry         (s3.py)
    |-- T3Geometry         (t3.py)
    |-- Nil3Geometry       (nil3.py)
    `-- Sol3Geometry       (sol3.py)
```

`UnifiedEngine` selects the concrete implementation from a `DOFConfig`.

---

## Modules

### `unified.py` - UnifiedEngine and DOFConfig

Primary runtime entry point.

**Key classes and functions:**

- `TopologyType`: enum with `S3`, `T3`, `NIL3`, `SOL3`.
- `FiberMode`: enum with `NONE`, `TWIST`, `MIXING`, `BOTH`.
- `CurvatureSource`: enum with `LC`, `EC`.
- `DOFConfig`: explicit configuration dataclass.
- `make_config(...)`: build a `DOFConfig` from topology/configuration fields.
- `make_engine(...)`: build the concrete topology engine.
- `normalize_topology_str(label)`: normalize labels to `S3`, `T3`, `Nil3`, or `Sol3`.

**Lorentzian LZ-native usage:**

Use `signature="lorentzian"` and `frame_convention="lz_native"` for Lorentzian LZ-native computations.

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

Abstract base class implementing topology-specific hooks for the engine pipeline.

**Subclass hooks:**

- `_build_radial_and_deformation_params(params)`: define scale/deformation symbols.
- `_build_structure_constants(params, C)`: fill frame structure constants.
- `_compute_volume(params)`: return the topology volume factor.

### `s3.py` - S3Geometry

`S3 x S1` engine based on the SU(2) Lie group.

- LZ reduced-sector specialization: `(a,b,c,u,v) = (4,4,4,0,0)`.
- LZ scalar `chi`: `4`.
- Closed isotropic reference topology.

### `t3.py` - T3Geometry

`T3 x S1` engine for the abelian flat three-torus.

- LZ reduced-sector specialization: `(a,b,c,u,v) = (0,0,0,0,0)`.
- LZ scalar `chi`: `0`.
- Flat reference topology.

### `nil3.py` - Nil3Geometry

`Nil3 x S1` engine for the Heisenberg nilmanifold.

- LZ reduced-sector specialization: `(a,b,c,u,v) = (0,0,-1,0,0)`.
- LZ scalar `chi`: `-1/12`.
- Supports a local isotropic-scale reduction.
- General Koszul formula is used because `Nil3` is not bi-invariant.

### `sol3.py` - Sol3Geometry

`Sol3 x S1` engine for the solvable Thurston geometry.

- LZ reduced-sector specialization: `(a,b,c,u,v) = (0,0,0,1,-1)`.
- LZ scalar `chi`: `-1/3`.
- Supports a local isotropic-scale reduction.
- Spatial structure constants are rigid under the diagonal frame scaling used here.

### `lz_invariants.py`

LZ-native topology invariant extraction.

**Key objects and functions:**

- `A, B, C, U, W`: five-parameter LZ-native symbols.
- `Q`: positive isotropic scale symbol.
- `CANONICAL_TOPOLOGIES`: `("S3", "T3", "Nil3", "Sol3")`.
- `engine_lz_structure_constants(topology)`: return LZ-native dense structure constants and the engine scale symbol.
- `extract_lz_five_parameter_specialization(topology)`: extract `(a,b,c,u,v)`.
- `topology_lz_parameter_table()`: return all canonical specializations.

This module is used by chi/scalar-curvature derivations and reduced-sector helpers.

### `lz_adapter.py`

Legacy-to-LZ index conversion helpers and verification utilities for earlier convention transitions.

### `bianchi_ix.py`

Auxiliary Bianchi IX / anisotropy helper used for comparison and guardrail calculations.  It is not part of the main `UnifiedEngine` topology dispatch.

---

## Topology Comparison

| Topology | LZ specialization `(a,b,c,u,v)` | `chi` | note |
|----------|---------------------------------|------:|----------------|
| `S3` | `(4,4,4,0,0)` | `4` | closed isotropic reference |
| `T3` | `(0,0,0,0,0)` | `0` | flat reference |
| `Nil3` | `(0,0,-1,0,0)` | `-1/12` | local isotropic-scale reduction |
| `Sol3` | `(0,0,0,1,-1)` | `-1/3` | local isotropic-scale reduction |

## Dependencies

- [engine](../engine/README.md): `BaseFrameEngine`, metric, connection, and pipeline infrastructure.
- [curvature](../curvature/README.md): curvature computation and spatial curvature helpers.
- [torsion](../torsion/README.md): torsion mode construction.
- [action](../action/README.md): action and effective-potential helpers.
