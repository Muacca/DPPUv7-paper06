# Utils Layer

-> [Japanese](README_ja.md)

Shared utility functions for tensor algebra, symbolic verification, logging, artifact I/O, and formatted output.

## Overview

The utils layer is intentionally small and dependency-light.  It provides:

- `tee_logger.py` for timestamped script logs.
- `epsilon.py` for metric-aware epsilon tensor contractions.
- `symbolic.py` for zero-proof helpers and tensor inventory checks.
- `io.py` for UTF-8 artifact I/O when CSV/JSON/text helpers are needed.

## Modules

### `epsilon.py`

Totally antisymmetric epsilon symbols and metric-aware raised/lowered epsilon tensors.

**Key functions:**

- `epsilon_symbol(a, b, c, d)`: 4D Levi-Civita symbol with lower-index orientation.
- `epsilon_4d(mu, nu, rho, sigma)`: 4D epsilon symbol.
- `epsilon_3d(i, j, k)`: 3D epsilon symbol.
- `epsilon_nd(indices)`: n-dimensional epsilon symbol.
- `metric_signature_sign(metric)`: sign of the metric determinant.
- `epsilon_tensor_down(metric, a, b, c, d)`: lowered epsilon tensor.
- `epsilon_tensor_up(metric_inv, a, b, c, d)`: raised epsilon tensor.

`epsilon_tensor_up` is useful for metric-aware epsilon contractions in symbolic diagnostics.

### `symbolic.py`

Symbolic simplification, zero-proving, witness search, and tensor inventory helpers.

**Key functions:**

- `prove_zero(expr, assumptions_dict=None, timeout_seconds=10)`: attempt symbolic zero proof.
- `generate_test_points(symbols_list, n_points=10)`: generate numerical test points.
- `find_nonzero_witness(expr, symbols_list, n_points=10, precision=50)`: search for numerical nonzero witnesses.
- `normalize_expression(expr)`: canonical simplification pass.
- `derivative_inventory(expr, symbol_pairs)`: report which symbols appear in an expression.
- `nonzero_component_count(array_expr, subs=None)`: count nonzero tensor components.

### `logger.py`

Step-oriented logger for the engine pipeline.

**Key classes:**

- `ComputationLogger`: log progress to file and console.
- `NullLogger`: suppress log output for tests or batch use.

`ComputationLogger` and `NullLogger` are re-exported from `dppu.engine` for backward compatibility.

### `tee_logger.py`

Tee-style stdout logger for scripts.

**Key functions/classes:**

- `setup_log(script_file, log_dir="logs")`: redirect stdout to both console and a timestamped log file.
- `teardown_log()`: restore stdout and print runtime metadata.
- `LogTee`: context-manager wrapper around `setup_log` / `teardown_log`.

Environment variables:

- `DPPU_LOG_DIR`: override the log directory.
- `DPPU_LOG_STDOUT=0`: write only to the log file and suppress console echo.

Callers choose the log directory through `log_dir` or `DPPU_LOG_DIR`.

### `io.py`

UTF-8 artifact I/O helpers.

**Key functions:**

- `now_utc_iso()`: UTC ISO-8601 timestamp.
- `read_text(path)` / `write_text(path, text)`: UTF-8 text I/O.
- `read_json(path)` / `write_json(path, payload)`: deterministic UTF-8 JSON I/O.
- `read_csv(path, missing_ok=False)` / `write_csv(path, fieldnames, rows)`: CSV helpers.
- `upsert_csv_rows(path, fieldnames, rows, key_fields)`: stable keyed CSV upsert.

### `printing.py`

Formatted console output helpers.

**Key functions:**

- `hline(char="=", width=70)`: horizontal separator.
- `print_header(title, width=70)`: section header.
- `print_sub(title)`: subsection label.

## Usage

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

## Dependencies

- SymPy for symbolic expressions.
- NumPy for numerical witness generation.

## Related Modules

- [curvature](../curvature/README.md): Hodge and Pontryagin contractions.
- [torsion](../torsion/README.md): torsion ansatz construction.
- [engine](../engine/README.md): pipeline logging and checkpoint integration.
