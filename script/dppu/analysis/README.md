# Analysis Layer

-> [Japanese](README_ja.md)

Symbolic classification and diagnostic-comparison helpers for DPPU-LZ quantities.

## Overview

The analysis layer provides utilities for classifying how symbolic quantities change under signature, metric, Hodge, index-convention, and torsion-ansatz choices.  It is a reporting and diagnostic layer: it does not decide physical stability, protection, or model validity.

## Modules

### `wick_classification.py`

Classification schema and diagnostic table helpers.

**Key constants:**

- `CLASSIFICATION_LABELS`: recognized behavior labels such as `INVARIANT`, `SIGN_FLIPPED`, `METRIC_DEPENDENT`, `HODGE_DEPENDENT`, and `TORSION_ANSATZ_DEPENDENT`.
- `CONFIDENCE_LABELS`: confidence labels such as `PROVEN_SYMBOLIC`, `CHECKED_BY_SCRIPT`, `DIAGNOSTIC_ONLY`, `HEURISTIC`, and `UNRESOLVED`.

**Key functions:**

- `classify_quantity_wick_behavior(quantity_name, evidence)`: return the classification row for a named quantity, optionally appending evidence.
- `build_phase1e_classification_table()`: return the bundled quantity-level classification table.
- `compare_pontryagin_diagnostic_cases(cases)`: return diagnostic comparison rows for recognized Pontryagin cases or unresolved placeholders for unknown cases.
- `render_classification_table_markdown(rows)`: render quantity classification rows as a Markdown table.
- `render_diagnostic_table_markdown(rows)`: render diagnostic comparison rows as a Markdown table.

The diagnostic case helpers are marked as diagnostic-only.  Their outputs are for algebraic comparison and should not be treated as physical models.

## Usage

```python
from dppu.analysis import (
    build_phase1e_classification_table,
    render_classification_table_markdown,
)

rows = build_phase1e_classification_table()
markdown = render_classification_table_markdown(rows)
```

```python
from dppu.analysis.wick_classification import compare_pontryagin_diagnostic_cases

rows = compare_pontryagin_diagnostic_cases([
    {"case_id": "LZ", "setup": "adopted Lorentzian metric-aware setup"},
])
```

## Dependencies

- SymPy for symbolic tensor expressions.
- [utils](../utils/README.md): epsilon symbols and metric-aware epsilon tensors.
- [curvature](../curvature/README.md): Pontryagin-density helper used by selected live diagnostics.

## Notes

- This layer classifies algebraic behavior only.
- `DIAGNOSTIC_ONLY` rows are comparison artifacts, not physical model definitions.
- Functions with historical names are kept for API compatibility.
