# Logging Conventions

This document defines the logging architecture and conventions for 
verification scripts under `scripts/`.

---

## 1. Two-Layer Architecture

The logging system has two independent layers with distinct responsibilities.

| Layer | API | Who uses it |
|-------|-----|-------------|
| **Engine layer** | `ComputationLogger` / `NullLogger` | `dppu/` library internals |
| **Script layer** | `setup_log` / `teardown_log` | Verification scripts under `scripts/` |

These two layers are completely independent.  
verification scripts use the script layer for user-facing output. When
they construct the engine directly, they rely on the engine default `NullLogger`
or pass `NullLogger()` explicitly; they do not use `ComputationLogger` directly.

Figure-rendering scripts are the intentional
exception: they write image files via `--output` and do not create execution
logs.

---

## 2. Engine Layer (`ComputationLogger` / `NullLogger`)

Located in `dppu/utils/logger.py`.

- `ComputationLogger`: writes pipeline step timings to a dedicated log file.
- `NullLogger`: no-op sink; discards all log calls.

**Rule**: Verification scripts rely on `NullLogger()` when constructing
`UnifiedEngine` or similar objects. Script output is handled entirely by the
script layer.

```python
from dppu.utils.logger import NullLogger

eng = UnifiedEngine(cfg, NullLogger())
```

---

## 3. Script Layer (`setup_log` / `teardown_log`)

Located in `dppu/utils/tee_logger.py`.

`setup_log` replaces `sys.stdout` with a `_Tee` that writes to both the
console and a log file simultaneously.  Every subsequent `print()` call is
captured automatically.  `teardown_log` restores `sys.stdout`, closes the log
file, and prints a summary line.

### 3.1 What `teardown_log` prints

```
(blank line)
Computation time: X.Xs
(blank line)
*Script: scripts/<subdir>/<name>.py
```

Scripts that already print their own timing or `*Script:` line should remove
them to avoid duplication.

### 3.2 Log file location

| Script type | `log_dir` argument |
|-------------|-------------------|
| argparse verification scripts | `args.output_dir` |
| non-argparse scripts | `"output"` (default) |

The log file is placed in the specified directory (created automatically) with
the name `<script_basename>_YYYYMMDD_HHMMSS.log`.

### 3.3 Environment variable overrides

| Variable | Effect |
|----------|--------|
| `DPPU_LOG_DIR` | Override `log_dir` with this path |
| `DPPU_LOG_STDOUT` | Set to `0` / `false` to suppress console echo (log file only) |

---

## 4. Usage Patterns

### Pattern A — argparse script with `main()`

```python
import argparse
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from dppu.utils.tee_logger import setup_log, teardown_log

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()

    setup_log(__file__, log_dir=args.output_dir)

    # ... computation using print() ...

    teardown_log()
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
```

### Pattern B — non-argparse script with `main()`

```python
if __name__ == "__main__":
    import os as _os, sys as _sys
    _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..'))
    from dppu.utils.tee_logger import setup_log, teardown_log
    setup_log(__file__, log_dir='output')
    main()
    teardown_log()
```

The `import` is intentionally placed inside the `if __name__` block so that
`sys.path` modification and `dppu` imports only occur when the script is run
directly (not when imported as a module).

### Pattern C — module-level script (no `main()`)

```python
import sys, os
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')  # encoding fix first

import time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from dppu.utils.tee_logger import setup_log, teardown_log
setup_log(__file__, log_dir='output')   # must come after encoding fix

# ... computation at module level using print() ...

teardown_log()   # at the very end of the file
```

**Important**: `setup_log` must be called **after** the Windows encoding fix
(`sys.stdout = io.TextIOWrapper(...)`) so that the `_Tee` wraps the
already-corrected UTF-8 stdout.

---

## 5. Output Conventions

- `print()` is the **sole** output method in scripts.  Do not use the stdlib
  `logging` module.
- Section separators use `"=" * 72` in Paper05 verification scripts.
- Step headers use `"--- Step N: ..."` or similar freeform text.
- `scripts/visualize/` may print a single completion line such as
  `[OK] wrote <path>` and does not need `setup_log`.

---

## 6. File References

| File | Role |
|------|------|
| `dppu/utils/tee_logger.py` | `setup_log`, `teardown_log`, `LogTee` |
| `dppu/utils/logger.py` | `ComputationLogger`, `NullLogger` |
