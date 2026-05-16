# Forms Layer

-> [Japanese](README_ja.md)

Exterior differential form algebra and Nieh-Yan form-engine helpers.

## Overview

The forms layer provides a lightweight symbolic exterior algebra representation plus a self-contained Nieh-Yan form engine.  It is separate from the component-tensor engine pipeline: this package works directly with forms represented as dictionaries of ordered basis-index tuples.

## Form Representation

A form is represented as:

```python
Form = dict[tuple[int, ...], object]
```

Examples:

- `{(): c}` is a scalar 0-form.
- `{(0,): 1}` is `e^0`.
- `{(0, 1): c}` is `c e^0 wedge e^1`.

Index tuples are stored in ascending order, and wedge-product signs are absorbed into the coefficient.

## Modules

### `exterior_algebra.py`

Pure exterior algebra operations.

**Key objects and functions:**

- `Form`: type alias for the dictionary representation.
- `basis(i)`: return the elementary 1-form `e^i`.
- `scalar_form(value)`: wrap a scalar as a 0-form.
- `clean_form(form)`: remove zero-coefficient entries after simplification.
- `add_forms(*forms)`: add forms.
- `scale_form(value, form)`: multiply a form by a scalar.
- `wedge(left, right)`: compute the exterior product.
- `form_to_str(form)`: human-readable representation.
- `coefficient(form, key=(0, 1, 2, 3))`: extract a coefficient.

### `nieh_yan.py`

Exterior-form engine for Nieh-Yan boundary forms and exact-form checks.

**Key class:**

- `NiehYanFormEngine`: stores a symbol set and computes torsion forms, Nieh-Yan 3-form `Q`, `dQ`, `T^a wedge T_a`, `e^a wedge e^b wedge R_ab`, and related audit data.

**Key methods:**

- `structure_constants()`: build frame structure constants for the engine's coframe convention.
- `exterior_derivative(form, C=None)`: compute exterior derivative.
- `torsion_tensor_low(mode)`: build lowered torsion tensor for `AX`, `VT`, or `MX`.
- `torsion_forms_low(mode)` / `torsion_forms_up(mode)`: build torsion 2-forms.
- `q_nieh_yan(mode)`: compute `Q = e^a wedge T_a`.
- `integrated_q_boundary(mode)`: compute the integrated boundary functional.
- `dq_reduced_density(mode)`: compute reduced density from `dQ`.
- `tt_form(mode)`: compute `T^a wedge T_a`.
- `curvature_forms_low(mode)`: compute EC curvature 2-forms.
- `e_e_r_form(mode)`: compute `e^a wedge e^b wedge R_ab`.
- `ny_form_audit(mode)`: return a dictionary comparing `dQ` with `TT - e e R`.

## Usage

```python
from sympy import symbols
from types import SimpleNamespace

from dppu.forms import NiehYanFormEngine, basis, form_to_str, wedge

e01 = wedge(basis(0), basis(1))
print(form_to_str(e01))

r, rd, rdd, N, Nd, eta, etad, V, Vd, V3, kappa, theta_NY = symbols(
    "r rd rdd N Nd eta etad V Vd V3 kappa theta_NY"
)
sy = SimpleNamespace(
    r=r, rd=rd, rdd=rdd, N=N, Nd=Nd,
    eta=eta, etad=etad, V=V, Vd=Vd,
    V3=V3, kappa=kappa, theta_NY=theta_NY,
)

engine = NiehYanFormEngine(sy)
audit = engine.ny_form_audit("MX")
```

## Relationship to `dppu.torsion.nieh_yan`

There are two Nieh-Yan-related modules by design:

- `dppu.torsion.nieh_yan` defines the `NyVariant` enum used to select `TT`, `REE`, or `FULL` component-density variants in the engine pipeline.
- `dppu.forms.nieh_yan` defines `NiehYanFormEngine`, an exterior-form engine for `Q`, `dQ`, and exact-form diagnostics.

The first module is a variant-selection API; the second is a form-calculation engine.

## Dependencies

- SymPy for symbolic expressions and simplification.
- [utils](../utils/README.md): epsilon symbols.
- [torsion](../torsion/README.md): torsion mode naming convention.
