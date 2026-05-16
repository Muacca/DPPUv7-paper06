"""Derive the paper06 auxiliary-shell reduced orbit atlas.

For each topology and mode, the script obtains F_branch from the same
engine-backed reduced-sector helper used by admissibility_classification.py,
solves the auxiliary equations, derives qdot^2 + F_on_shell = 0, and classifies
the real orbit sheets inside the paper06 bundle.

Outputs:
  data/orbit_atlas_YYYYMMDD_HHMMSS.log
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
SCRIPT_ROOT = HERE.parents[1]     # .../script/
PROJECT_ROOT = HERE.parents[2]    # .../20_DPPUv7-paper06/
DATA_DIR = PROJECT_ROOT / "data"

if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from dppu.utils.tee_logger import setup_log, teardown_log

from dppu.action.reduced_sector import (
    MODES,
    QDOT,
    TOPOLOGIES,
    chi_for_topology,
    derive_static_branch_function,
    solve_auxiliary_shell,
)


NO_REAL_AUX_SHELL_ORBIT = "NO_REAL_AUX_SHELL_ORBIT"
STATIC_OR_DEGENERATE = "STATIC_OR_DEGENERATE"
MONOTONIC_EXPANSION = "MONOTONIC_EXPANSION"
SINGULAR_APPROACH = "SINGULAR_APPROACH"
BOUNCE_LIKE = "BOUNCE_LIKE"
RECOLLAPSE_LIKE = "RECOLLAPSE_LIKE"

Q0 = 1.0
T_SPAN = 1.0
N_STEPS = 200
Q_FLOOR = 1e-8

ROW_FIELDS = [
    "topology", "mode", "chi_from_lc_curvature", "F_on_shell",
    "orbit_equation", "sheet", "qdot0_analytic", "orbit_class",
    "constraint_residual_max", "constraint_check", "overall_status",
]


def format_expr(expr) -> str:
    if expr is None:
        return "n/a"
    if isinstance(expr, list):
        return ",".join(str(item) for item in expr) if expr else "none"
    if isinstance(expr, dict):
        return ",".join(f"{key}={value}" for key, value in expr.items()) if expr else "trivial"
    return str(sp.factor(sp.cancel(expr))) if isinstance(expr, sp.Basic) else str(expr)


def integrate_linear_orbit(chi_value: float, qdot0: float) -> float:
    """Integrate q''=0 and return max |qdot^2 + chi|."""
    dt = T_SPAN / N_STEPS
    q = Q0
    qdot = qdot0
    residuals = [abs(qdot**2 + chi_value)]
    for _ in range(N_STEPS):
        q += dt * qdot
        residuals.append(abs(qdot**2 + chi_value))
        if q <= Q_FLOOR:
            break
    return float(max(residuals))


def rows_for_branch(topology: str, mode: str) -> list[dict[str, str]]:
    chi_value = chi_for_topology(topology)
    branch_f = derive_static_branch_function(topology, mode)
    shell = solve_auxiliary_shell(branch_f)
    f_on_shell = sp.factor(sp.cancel(sp.simplify(shell["on_shell"])))
    orbit_equation = sp.factor(sp.cancel(QDOT**2 + f_on_shell))
    sign = sp.signsimp(f_on_shell)

    base = {
        "topology": topology,
        "mode": mode,
        "chi_from_lc_curvature": format_expr(chi_value),
        "F_on_shell": format_expr(f_on_shell),
        "orbit_equation": format_expr(orbit_equation),
    }

    if sign > 0:
        return [
            {
                **base,
                "sheet": "none",
                "qdot0_analytic": "no_real_solution",
                "orbit_class": NO_REAL_AUX_SHELL_ORBIT,
                "constraint_residual_max": "n/a",
                "constraint_check": "PASS",
                "overall_status": "PASS",
            }
        ]

    if sign == 0:
        residual = sp.factor(sp.cancel(orbit_equation.subs(QDOT, 0)))
        status = "PASS" if residual == 0 else "FAIL"
        return [
            {
                **base,
                "sheet": "zero",
                "qdot0_analytic": "0",
                "orbit_class": STATIC_OR_DEGENERATE,
                "constraint_residual_max": format_expr(residual),
                "constraint_check": status,
                "overall_status": status,
            }
        ]

    qdot_mag = sp.sqrt(-f_on_shell)
    rows = []
    for sheet, qdot_exact, orbit_class in [
        ("expanding", qdot_mag, MONOTONIC_EXPANSION),
        ("collapsing", -qdot_mag, SINGULAR_APPROACH),
    ]:
        residual_exact = sp.factor(sp.cancel(orbit_equation.subs(QDOT, qdot_exact)))
        residual_numeric = integrate_linear_orbit(float(f_on_shell), float(qdot_exact))
        status = "PASS" if residual_exact == 0 and residual_numeric < 1e-12 else "FAIL"
        rows.append(
            {
                **base,
                "sheet": sheet,
                "qdot0_analytic": format_expr(qdot_exact),
                "orbit_class": orbit_class,
                "constraint_residual_max": f"{residual_numeric:.6e}",
                "constraint_check": status,
                "overall_status": status,
            }
        )
    return rows


def print_rows(rows: list[dict[str, str]]) -> None:
    print("row_details_begin")
    for row in rows:
        ordered = {field: row[field] for field in ROW_FIELDS}
        print("  row_detail " + json.dumps(ordered, sort_keys=False, ensure_ascii=True))
    print("row_details_end")


def main() -> int:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    setup_log(__file__, log_dir=str(DATA_DIR))
    try:
        print("Paper06 reduced orbit atlas derivation")
        print("=" * 65)
        rows = []
        for topology in TOPOLOGIES:
            for mode in MODES:
                branch_rows = rows_for_branch(topology, mode)
                rows.extend(branch_rows)
                for row in branch_rows:
                    print(
                        f"  {topology}x{mode}: "
                        f"sheet={row['sheet']} "
                        f"class={row['orbit_class']} "
                        f"constraint={row['constraint_check']} "
                        f"status={row['overall_status']}"
                    )

        print_rows(rows)

        class_counts: dict[str, int] = {}
        for row in rows:
            label = row["orbit_class"]
            class_counts[label] = class_counts.get(label, 0) + 1
        fail_count = sum(1 for row in rows if row["overall_status"] != "PASS")
        bounce_count = class_counts.get(BOUNCE_LIKE, 0) + class_counts.get(RECOLLAPSE_LIKE, 0)

        print()
        print(f"total_orbit_entries={len(rows)}")
        print(f"class_counts={class_counts}")
        print(f"bounce_or_recollapse_count={bounce_count}")
        print(f"fail_count={fail_count}")
        print(f"overall_verdict={'PASS' if fail_count == 0 else 'FAIL'}")
        return 0 if fail_count == 0 else 1
    finally:
        teardown_log()


if __name__ == "__main__":
    raise SystemExit(main())
