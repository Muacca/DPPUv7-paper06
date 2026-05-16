# Engine Layer

-> [Japanese](README_ja.md)

Computation pipeline infrastructure and low-level frame-geometry helpers.

## Overview

The engine layer provides the base pipeline used by topology engines, plus utilities for metric handling, Levi-Civita connection construction, contortion, Einstein-Cartan connection construction, spin-2 helpers, and checkpointing.

`ComputationLogger` and `NullLogger` live in [utils](../utils/README.md) and are re-exported from `dppu.engine` for backward compatibility.

## Modules

### `pipeline.py`

Base class for topology engines.

**BaseFrameEngine step order:**

1. `E4.1` - parameter setup.
2. `E4.2` - metric, frame, volume, and structure constants.
3. `E4.3` - Levi-Civita connection.
4. `E4.3a` - LC Riemann/Ricci path.
5. `E4.3b` - LC Weyl tensor path.
6. `E4.4` - torsion ansatz.
7. `E4.5` - contortion tensor.
8. `E4.6` - Einstein-Cartan connection.
9. `E4.7` - EC Riemann tensor.
10. `E4.8` - EC Ricci scalar.
11. `E4.8b` - EC Weyl scalar.
12. `E4.9` - torsion scalar.
13. `E4.10` - Nieh-Yan density.
14. `E4.11` - Lagrangian construction.
15. `E4.12` - angular/topology volume integration.
16. `E4.13` - effective-potential extraction.
17. `E4.15` - summary.

The `run(start_step="E4.1")` method runs from the requested start step to the end.  Call individual `step_E4_*` methods directly when only a partial pipeline is needed.

### `metric.py`

Frame metric utilities.

### `levi_civita.py`

Levi-Civita connection from frame structure constants, using the Koszul formula.

### `contortion.py`

Contortion tensor construction from torsion.

### `ec_connection.py`

Einstein-Cartan connection construction from the LC connection and contortion.

### `spin2.py`

Spin-2/off-diagonal deformation helpers used by higher-order diagnostics.

### `checkpoint.py`

Save and restore computation state.

```python
from dppu.engine import CheckpointManager

ckpt = CheckpointManager("checkpoints/", enabled=True)
ckpt.save("E4.7", {"R": riemann_tensor})
data = ckpt.load("E4.7")
```

## Usage

```python
from dppu.topology import FiberMode, TopologyType, make_engine
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

engine = make_engine(
    TopologyType.S3,
    torsion_mode=Mode.MX,
    ny_variant=NyVariant.FULL,
    signature="lorentzian",
    frame_convention="lz_native",
    fiber_mode=FiberMode.NONE,
    enable_squash=False,
)
engine.run()
```

Partial execution can be done by calling explicit step methods:

```python
engine.step_E4_1_setup()
engine.step_E4_2_metric_and_frame()
engine.step_E4_3_christoffel_frame()
engine.step_E4_4_torsion_ansatz_frame()
engine.step_E4_5_contortion_frame()
engine.step_E4_6_ec_connection_frame()
engine.step_E4_7_riemann_tensor_frame()
```

## Dependencies

- SymPy for symbolic tensor expressions.
- [utils](../utils/README.md): logging, epsilon tensors, and symbolic helpers.

## Related Modules

- [topology](../topology/README.md): concrete topology engines.
- [torsion](../torsion/README.md): torsion ansatz construction.
- [curvature](../curvature/README.md): Riemann, Ricci, Weyl, and Pontryagin routines.
- [action](../action/README.md): Lagrangian and effective-potential routines.
