# Curvature Layer

-> [Japanese](README_ja.md)

Module group for curvature tensor computation, Pontryagin diagnostics, Hodge operations, Weyl tensors, and spatial-curvature helpers.

## Overview

The curvature layer provides symbolic tensor routines used by the EC+NY engine and downstream callers.  The layer includes:

- `spatial_lie.py`: derives the scale-stripped 3D scalar curvature and `chi`.
- `pontryagin_lz.py`: computes the Lorentzian form-Hodge Pontryagin density used for `P_form` checks.

These modules keep form-Hodge Pontryagin densities separate from any script-level diagnostic contractions.

## Modules

### `riemann.py`

Riemann tensor computation in a non-coordinate frame, including strict antisymmetry checks and first-index lowering.

**Key classes/functions:**

- `RiemannAntisymmetryError`: exception raised by strict verification.
- `verify_antisymmetry_strict(R, dim, ...)`: verify Riemann antisymmetry identities.
- `compute_riemann_tensor(gamma, structure_constants, dim, ...)`: compute `R^a_{bcd}`.
- `lower_first_index(Riemann, metric, dim)`: lower the first index to obtain `R_abcd`.

### `ricci.py`

Ricci tensor, Ricci scalar, and decomposition helpers.

**Key functions:**

- `compute_ricci_tensor(Riemann, dim)`: Ricci contraction.
- `compute_ricci_scalar(Riemann, dim, logger=None, metric_inv=None)`: Ricci scalar from Riemann data.
- `compute_ricci_scalar_from_tensor(Ricci, metric_inv, dim)`: scalar from a Ricci tensor.
- `decompose_ricci_tensor(Ricci, metric, dim)`: trace/traceless decomposition.

### `hodge.py`

Hodge dual and 2-form block helpers.

**Key functions:**

- `hodge_dual_2form(F_down, metric, metric_inv=None, signature="lorentzian")`: Hodge dual for a 2-form.
- `basis_2form(a, b)`: basis 2-form.
- `compute_hodge_dual(R_tensor)`: legacy Riemann Hodge dual helper.
- `classify_lz_2form_block(c, d)`: LZ-native block classification.
- `hodge_swaps_blocks()`: block-swap sanity check.

### `pontryagin.py`

Legacy numerical Pontryagin helper.

**Key functions:**

- `compute_P_from_riemann(R_arr)`: compute `P=<R,*R>` for numerical arrays.
- `get_riemann_numeric(...)`: extract numerical Riemann data from a lambdified engine result.

### `pontryagin_lz.py`

Lorentzian/LZ-native Pontryagin density helpers.

**Key functions:**

- `pontryagin_density_lorentzian(R_abcd, metric, simplify_result=True)`: compute the form-Hodge `P_form=<R,*R>` density.
- `pontryagin_density_lorentzian_by_epsilon(...)`: epsilon-contraction cross-check route.
- `classify_pontryagin_expression(expr, zero_tol=None)`: classify exact-zero/nonzero symbolic results.
- `classify_pontryagin_expression_deep(expr, zero_tol=None)`: deeper simplification and classification route.

This module is intended for form-Hodge `P_form=<R,*R>` computations and exact-zero classification.  It is separate from any custom internal-pair diagnostic contractions implemented by individual scripts.

### `spatial_lie.py`

Scale-stripped spatial Lie-frame curvature helpers.

**Key functions:**

- `spatial_structure_constants_3d(params=None)`: build 3D structure constants for the five-parameter LZ family.
- `scalar_curvature_q2_from_structure_constants(structure_constants)`: compute `R^(3) q^2`.
- `chi_from_lz_parameters(params)`: compute `chi = R^(3) q^2 / 6`.

### `sd_extension.py`

Self-duality (SD/ASD) mixin for engines.

**Key classes:**

- `SDExtensionMixin`: attaches SD functionality to topology engines.

### `sd_diagnostics.py`

Self-duality and Pontryagin diagnostic evaluation.

**Key classes/functions:**

- `CurvatureSDDiagnostics`: evaluate SD/ASD residuals and related curvature diagnostics.

### `weyl.py`

Weyl tensor and Weyl scalar computation.

**Key functions:**

- `compute_weyl_tensor(R_abcd, Ricci, R_scalar, metric, dim)`: compute `C_abcd`.
- `compute_weyl_scalar(C_abcd, metric)`: compute `C^2`.

## Notes

- `P_form` is the form-Hodge density computed through `pontryagin_lz.py`.
- Internal-pair diagnostic contractions are script-level quantities and should not be conflated with the true form-Hodge Pontryagin density.
- `chi` is derived from the LC/spatial-curvature sector through `spatial_lie.py`.

## Dependencies

- [engine](../engine/README.md): LC/EC connection data and frame pipeline.
- [topology](../topology/README.md): LZ-native structure constants.
- [utils](../utils/README.md): epsilon tensors and symbolic utilities.
