# Action Layer

-> [English](README.md)

Einstein-Cartan action、effective potential、および reduced-sector branch extraction のためのモジュール群です。

## 概要

action layer は topology engine と downstream callers が利用する symbolic action utility を提供します。reduced-sector の入口は `reduced_sector.py` です。このモジュールは static branch function `F_branch` を抽出し、topology data から `chi` を導出し、`EH/AX/VT/MX` branches の auxiliary torsion shell を解きます。

## モジュール

### `lagrangian.py`

Einstein-Cartan + Nieh-Yan Lagrangian の構成。

**Action:**

```text
S = integral d^4x sqrt(|g|) L
L = R/(2 kappa^2) + theta_NY * N + alpha * C^2
```

**主な関数:**

- `compute_lagrangian(ricci_scalar, nieh_yan_density, kappa, theta_NY, weyl_scalar=0, alpha=0, logger=None)`: EC+NY(+Weyl) Lagrangian density を構成します。
- `compute_action(L_density, volume_factor)`: 与えられた volume factor を掛けて action expression を構成します。

### `potential.py`

Effective potential の抽出と分解。

**主な関数:**

- `subs_zero_modes(expr, params)`: inactive symbols を zero に置換します。
- `compute_effective_potential(action, logger=None)`: `V_eff = -action` を計算します。
- `get_potential_function(potential, param_symbols)`: numerical function を構成します。
- `decompose_potential(potential, r)`: symbolic expression を有用な成分に分解します。

### `ec_action.py`

Engine pipeline で使う Einstein-Cartan effective-potential helper。

**主な関数:**

- `compute_c2_ec(data)`: EC curvature/torsion coefficient package を計算します。
- `build_veff_ec(topology_type, torsion_mode=None, ny_variant=None, ...)`: symbolic EC effective potential を構成します。
- `build_veff_ec_func(V_eff, symbols)`: effective-potential expression を lambdify します。

### `reduced_sector.py`

Reduced-sector extraction helper。

topology dependence を少数の symbolic branch functions と auxiliary equations に局在させるための層です。

**主な object / 関数:**

- `TOPOLOGIES`: canonical topology list、`("S3", "T3", "Nil3", "Sol3")`。
- `MODES`: reduced branch list、`("EH", "AX", "VT", "MX")`。
- `StaticBranchResult`: engine-derived branch result を保持する dataclass。
- `chi_for_topology(topology)`: LZ-native topology specialization から `chi` を導出します。
- `derive_static_branch_result(topology, mode)`: branch metadata 一式を返します。
- `derive_static_branch_function(topology, mode)`: `F_branch` のみを返します。
- `auxiliary_variables(expr)`: active auxiliary torsion variables を抽出します。
- `solve_auxiliary_shell(expr)`: auxiliary equations を解き、Hessian / on-shell data を返します。

呼び出し側では、これらの helper により次の reduced-sector identity を確認できます。

```text
F_branch | auxiliary shell = chi
qdot^2 + chi = 0
```

## 使用例

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

## 依存関係

- [curvature](../curvature/README_ja.md): spatial curvature と `chi` derivation。
- [topology](../topology/README_ja.md): topology engines と LZ-native structure constants。
- [torsion](../torsion/README_ja.md): torsion mode definitions と ansatz construction。
