# Engine Layer

-> [English](README.md)

Computation pipeline infrastructure と low-level frame-geometry helpers を提供します。

## 概要

engine layer は topology engines が利用する base pipeline と、metric handling、Levi-Civita connection construction、contortion、Einstein-Cartan connection construction、spin-2 helpers、checkpointing を提供します。

`ComputationLogger` と `NullLogger` は [utils](../utils/README_ja.md) にあり、backward compatibility のため `dppu.engine` からも re-export されています。

## モジュール

### `pipeline.py`

Topology engines の base class。

**BaseFrameEngine step order:**

1. `E4.1` - parameter setup。
2. `E4.2` - metric, frame, volume, structure constants。
3. `E4.3` - Levi-Civita connection。
4. `E4.3a` - LC Riemann/Ricci path。
5. `E4.3b` - LC Weyl tensor path。
6. `E4.4` - torsion ansatz。
7. `E4.5` - contortion tensor。
8. `E4.6` - Einstein-Cartan connection。
9. `E4.7` - EC Riemann tensor。
10. `E4.8` - EC Ricci scalar。
11. `E4.8b` - EC Weyl scalar。
12. `E4.9` - torsion scalar。
13. `E4.10` - Nieh-Yan density。
14. `E4.11` - Lagrangian construction。
15. `E4.12` - angular/topology volume integration。
16. `E4.13` - effective-potential extraction。
17. `E4.15` - summary。

`run(start_step="E4.1")` は指定 step から最後まで実行します。partial pipeline だけが必要な場合は、個別の `step_E4_*` methods を直接呼びます。

### `metric.py`

Frame metric utilities。

### `levi_civita.py`

Frame structure constants から Koszul formula によって Levi-Civita connection を構成します。

### `contortion.py`

Torsion から contortion tensor を構成します。

### `ec_connection.py`

LC connection と contortion から Einstein-Cartan connection を構成します。

### `spin2.py`

Higher-order diagnostics で使う spin-2/off-diagonal deformation helpers。

### `checkpoint.py`

Computation state の保存と復元。

```python
from dppu.engine import CheckpointManager

ckpt = CheckpointManager("checkpoints/", enabled=True)
ckpt.save("E4.7", {"R": riemann_tensor})
data = ckpt.load("E4.7")
```

## 使用例

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

Partial execution は explicit step methods の呼び出しで行います。

```python
engine.step_E4_1_setup()
engine.step_E4_2_metric_and_frame()
engine.step_E4_3_christoffel_frame()
engine.step_E4_4_torsion_ansatz_frame()
engine.step_E4_5_contortion_frame()
engine.step_E4_6_ec_connection_frame()
engine.step_E4_7_riemann_tensor_frame()
```

## 依存関係

- SymPy: symbolic tensor expressions。
- [utils](../utils/README_ja.md): logging, epsilon tensors, symbolic helpers。

## 関連モジュール

- [topology](../topology/README_ja.md): concrete topology engines。
- [torsion](../torsion/README_ja.md): torsion ansatz construction。
- [curvature](../curvature/README_ja.md): Riemann, Ricci, Weyl, Pontryagin routines。
- [action](../action/README_ja.md): Lagrangian and effective-potential routines。
