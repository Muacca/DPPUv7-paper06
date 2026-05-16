# Action Layer

-> [Japanese](README_ja.md)

Module group for Einstein-Cartan action, effective potential, and reduced-sector branch extraction.

## Overview

The action layer provides symbolic action utilities used by the topology engines and downstream callers.  The reduced-sector entry point is `reduced_sector.py`, which extracts the static branch function `F_branch`, derives `chi` from topology data, and solves the auxiliary torsion shell for the `EH/AX/VT/MX` branches.

## Modules

### `lagrangian.py`

Einstein-Cartan + Nieh-Yan Lagrangian construction.

**Action:**

```text
S = integral d^4x sqrt(|g|) L
L = R/(2 kappa^2) + theta_NY * N + alpha * C^2
```

**Key functions:**

- `compute_lagrangian(ricci_scalar, nieh_yan_density, kappa, theta_NY, weyl_scalar=0, alpha=0, logger=None)`: construct the EC+NY(+Weyl) Lagrangian density.
- `compute_action(L_density, volume_factor)`: multiply by the supplied volume factor.

### `potential.py`

Effective potential extraction and decomposition.

**Key functions:**

- `subs_zero_modes(expr, params)`: substitute inactive symbols by zero.
- `compute_effective_potential(action, logger=None)`: compute `V_eff = -action`.
- `get_potential_function(potential, param_symbols)`: build a numerical function.
- `decompose_potential(potential, r)`: separate useful symbolic pieces.

### `ec_action.py`

Einstein-Cartan effective-potential helpers used by the engine pipeline.

**Key functions:**

- `compute_c2_ec(data)`: compute the EC curvature/torsion coefficient package.
- `build_veff_ec(topology_type, torsion_mode=None, ny_variant=None, ...)`: build the symbolic EC effective potential.
- `build_veff_ec_func(V_eff, symbols)`: lambdify an effective-potential expression.

### `reduced_sector.py`

Reduced-sector extraction helper.

This module keeps topology dependence in a small set of symbolic branch functions and auxiliary equations.

**Key objects and functions:**

- `TOPOLOGIES`: canonical topology list, `("S3", "T3", "Nil3", "Sol3")`.
- `MODES`: reduced branch list, `("EH", "AX", "VT", "MX")`.
- `StaticBranchResult`: dataclass containing the engine-derived branch result.
- `chi_for_topology(topology)`: derive `chi` from the LZ-native topology specialization.
- `derive_static_branch_result(topology, mode)`: return full branch metadata.
- `derive_static_branch_function(topology, mode)`: return only `F_branch`.
- `auxiliary_variables(expr)`: identify active auxiliary torsion variables.
- `solve_auxiliary_shell(expr)`: solve auxiliary equations and return Hessian/on-shell data.

Callers can use these helpers to verify the reduced-sector identity

```text
F_branch | auxiliary shell = chi
qdot^2 + chi = 0
```

## Usage

```python
import sympy as sp

from dppu.action.reduced_sector import (
    TOPOLOGIES,
    MODES,
    chi_for_topology,
    derive_static_branch_function,
    solve_auxiliary_shell,
)

for topology in TOPOLOGIES:
    for mode in MODES:
        branch_f = derive_static_branch_function(topology, mode)
        shell = solve_auxiliary_shell(branch_f)
        assert sp.simplify(shell["on_shell"] - chi_for_topology(topology)) == 0
```

## Dependencies

- [curvature](../curvature/README.md): spatial curvature and `chi` derivation.
- [topology](../topology/README.md): topology engines and LZ-native structure constants.
- [torsion](../torsion/README.md): torsion mode definitions and ansatz construction.
