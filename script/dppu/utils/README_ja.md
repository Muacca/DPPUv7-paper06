# Utils Layer

-> [English](README.md)

Tensor algebra、symbolic verification、logging、artifact I/O、formatted output のための共通 utility 群です。

## 概要

utils layer は小さく、依存を軽く保つための共通層です。次の utility を提供します。

- `tee_logger.py`: timestamped script logs。
- `epsilon.py`: metric-aware epsilon tensor contractions。
- `symbolic.py`: zero-proof helpers と tensor inventory checks。
- `io.py`: CSV/JSON/text が必要な場合の UTF-8 artifact I/O。

## モジュール

### `epsilon.py`

Totally antisymmetric epsilon symbols と metric-aware raised/lowered epsilon tensors。

**主な関数:**

- `epsilon_symbol(a, b, c, d)`: lower-index orientation を持つ 4D Levi-Civita symbol。
- `epsilon_4d(mu, nu, rho, sigma)`: 4D epsilon symbol。
- `epsilon_3d(i, j, k)`: 3D epsilon symbol。
- `epsilon_nd(indices)`: n-dimensional epsilon symbol。
- `metric_signature_sign(metric)`: metric determinant の sign。
- `epsilon_tensor_down(metric, a, b, c, d)`: lowered epsilon tensor。
- `epsilon_tensor_up(metric_inv, a, b, c, d)`: raised epsilon tensor。

`epsilon_tensor_up` は metric-aware epsilon contraction を含む symbolic diagnostics で利用できます。

### `symbolic.py`

Symbolic simplification、zero-proving、witness search、tensor inventory helper。

**主な関数:**

- `prove_zero(expr, assumptions_dict=None, timeout_seconds=10)`: symbolic zero proof を試みます。
- `generate_test_points(symbols_list, n_points=10)`: numerical test points を生成します。
- `find_nonzero_witness(expr, symbols_list, n_points=10, precision=50)`: numerical nonzero witness を探索します。
- `normalize_expression(expr)`: canonical simplification pass。
- `derivative_inventory(expr, symbol_pairs)`: expression 内に現れる symbols を報告します。
- `nonzero_component_count(array_expr, subs=None)`: tensor の nonzero components を数えます。

### `logger.py`

Engine pipeline 用の step-oriented logger。

**主なクラス:**

- `ComputationLogger`: file と console に progress を記録します。
- `NullLogger`: tests や batch use のために log output を抑制します。

`ComputationLogger` と `NullLogger` は backward compatibility のため `dppu.engine` からも re-export されます。

### `tee_logger.py`

Script 用の tee-style stdout logger。

**主な関数 / クラス:**

- `setup_log(script_file, log_dir="logs")`: stdout を console と timestamped log file の両方へ redirect します。
- `teardown_log()`: stdout を復元し、runtime metadata を出力します。
- `LogTee`: `setup_log` / `teardown_log` の context-manager wrapper。

Environment variables:

- `DPPU_LOG_DIR`: log directory を上書きします。
- `DPPU_LOG_STDOUT=0`: console echo を抑制し、log file のみに出力します。

呼び出し側は `log_dir` または `DPPU_LOG_DIR` により log directory を指定できます。

### `io.py`

UTF-8 artifact I/O helpers。

**主な関数:**

- `now_utc_iso()`: UTC ISO-8601 timestamp。
- `read_text(path)` / `write_text(path, text)`: UTF-8 text I/O。
- `read_json(path)` / `write_json(path, payload)`: deterministic UTF-8 JSON I/O。
- `read_csv(path, missing_ok=False)` / `write_csv(path, fieldnames, rows)`: CSV helpers。
- `upsert_csv_rows(path, fieldnames, rows, key_fields)`: stable keyed CSV upsert。

### `printing.py`

Formatted console output helpers。

**主な関数:**

- `hline(char="=", width=70)`: horizontal separator。
- `print_header(title, width=70)`: section header。
- `print_sub(title)`: subsection label。

## 使用例

```python
from dppu.utils.tee_logger import setup_log, teardown_log

setup_log(__file__, log_dir="../data")
try:
    print("verification output")
finally:
    teardown_log()
```

```python
from dppu.utils.epsilon import epsilon_tensor_up

eps = epsilon_tensor_up(metric_inv, 0, 1, 2, 3)
```

## 依存関係

- SymPy: symbolic expressions。
- NumPy: numerical witness generation。

## 関連モジュール

- [curvature](../curvature/README_ja.md): Hodge and Pontryagin contractions。
- [torsion](../torsion/README_ja.md): torsion ansatz construction。
- [engine](../engine/README_ja.md): pipeline logging and checkpoint integration。
