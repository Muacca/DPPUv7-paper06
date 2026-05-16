# Analysis Layer

-> [English](README.md)

DPPU-LZ quantities の symbolic classification と diagnostic-comparison helpers を提供します。

## 概要

analysis layer は、signature、metric、Hodge、index convention、torsion ansatz の選択によって symbolic quantity がどのように変化するかを分類する utility を提供します。この layer は reporting / diagnostic layer であり、physical stability、protection、model validity を判定する層ではありません。

## モジュール

### `wick_classification.py`

Classification schema と diagnostic table helpers。

**主な定数:**

- `CLASSIFICATION_LABELS`: `INVARIANT`, `SIGN_FLIPPED`, `METRIC_DEPENDENT`, `HODGE_DEPENDENT`, `TORSION_ANSATZ_DEPENDENT` などの behavior labels。
- `CONFIDENCE_LABELS`: `PROVEN_SYMBOLIC`, `CHECKED_BY_SCRIPT`, `DIAGNOSTIC_ONLY`, `HEURISTIC`, `UNRESOLVED` などの confidence labels。

**主な関数:**

- `classify_quantity_wick_behavior(quantity_name, evidence)`: named quantity に対する classification row を返し、必要に応じて evidence を追記します。
- `build_phase1e_classification_table()`: bundled quantity-level classification table を返します。
- `compare_pontryagin_diagnostic_cases(cases)`: 認識済み Pontryagin cases の diagnostic comparison rows、または unknown cases の unresolved placeholders を返します。
- `render_classification_table_markdown(rows)`: quantity classification rows を Markdown table として出力します。
- `render_diagnostic_table_markdown(rows)`: diagnostic comparison rows を Markdown table として出力します。

diagnostic case helpers は diagnostic-only です。出力は algebraic comparison のためのものであり、physical model として扱わないでください。

## 使用例

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

## 依存関係

- SymPy: symbolic tensor expressions。
- [utils](../utils/README_ja.md): epsilon symbols と metric-aware epsilon tensors。
- [curvature](../curvature/README_ja.md): selected live diagnostics で使う Pontryagin-density helper。

## 注記

- この layer は algebraic behavior の分類のみを扱います。
- `DIAGNOSTIC_ONLY` rows は comparison artifacts であり、physical model definitions ではありません。
- historical names を持つ関数は API compatibility のため維持されています。
