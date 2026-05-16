"""Derive the paper06 admissibility classification.

The branch functions are obtained from the bundled EC+NY engine for AX/VT/MX
and from the Appendix A spatial curvature for EH.  The script then performs
the lapse and auxiliary algebra symbolically inside the paper06 bundle.

Outputs:
  data/admissibility_classification_YYYYMMDD_HHMMSS.log
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

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
    N,
    QDOT,
    TOPOLOGIES,
    chi_for_topology,
    derive_static_branch_function,
    solve_auxiliary_shell,
)


ROW_FIELDS = [
    "topology", "mode", "chi_from_lc_curvature", "F_branch_engine_derived",
    "hamiltonian_constraint", "aux_vars", "aux_equations",
    "shell_solution", "F_on_shell", "F_on_shell_minus_chi",
    "hessian_det", "hessian_status", "admissibility_label",
    "label_origin", "overall_status",
]


def format_expr(expr) -> str:
    if expr is None:
        return "n/a"
    if isinstance(expr, list):
        return ",".join(str(item) for item in expr) if expr else "none"
    if isinstance(expr, dict):
        return ",".join(f"{key}={value}" for key, value in expr.items()) if expr else "trivial"
    return str(sp.factor(sp.cancel(expr))) if isinstance(expr, sp.Basic) else str(expr)


def has_mixed_auxiliary_coupling(shell_data: dict[str, object]) -> bool:
    """Return whether the auxiliary Hessian contains eta/V mixing."""
    aux_vars = shell_data["aux_vars"]
    hessian = shell_data["hessian"]
    if hessian is None or len(aux_vars) < 2:
        return False
    for row_index in range(len(aux_vars)):
        for col_index in range(len(aux_vars)):
            if row_index != col_index and sp.simplify(hessian[row_index, col_index]) != 0:
                return True
    return False


def classify_admissibility(shell_data: dict[str, object], on_shell_minus_chi: sp.Expr) -> tuple[str, str]:
    """Classify paper06 branch admissibility from computed auxiliary data."""
    if on_shell_minus_chi != 0:
        return "L_INADMISSIBLE", "F_ON_SHELL_NOT_CHI"
    if shell_data["solution_status"] not in {"NO_AUXILIARY_VARIABLES", "UNIQUE_NONDEGENERATE"}:
        return "L_INADMISSIBLE", str(shell_data["solution_status"])
    if has_mixed_auxiliary_coupling(shell_data):
        return "L_CONDITIONALLY_ADMISSIBLE", "MIXED_AUXILIARY_HESSIAN"
    return "L_ADMISSIBLE", "AUXILIARY_CLOSURE"


def derive_row(topology: str, mode: str) -> dict[str, str]:
    chi_value = chi_for_topology(topology)
    branch_f = derive_static_branch_function(topology, mode)
    shell = solve_auxiliary_shell(branch_f)
    f_on_shell = shell["on_shell"]
    f_minus_chi = sp.factor(sp.cancel(sp.simplify(f_on_shell - chi_value)))
    label, label_origin = classify_admissibility(shell, f_minus_chi)
    hamiltonian_constraint = sp.factor(sp.cancel(QDOT**2 / N**2 + branch_f))
    status = "PASS" if label != "L_INADMISSIBLE" else "FAIL"

    return {
        "topology": topology,
        "mode": mode,
        "chi_from_lc_curvature": format_expr(chi_value),
        "F_branch_engine_derived": format_expr(branch_f),
        "hamiltonian_constraint": format_expr(hamiltonian_constraint),
        "aux_vars": format_expr(shell["aux_vars"]),
        "aux_equations": format_expr(shell["equations"]),
        "shell_solution": format_expr(shell["solution"]),
        "F_on_shell": format_expr(f_on_shell),
        "F_on_shell_minus_chi": format_expr(f_minus_chi),
        "hessian_det": format_expr(shell["hessian_det"]),
        "hessian_status": str(shell["hessian_status"]),
        "admissibility_label": label,
        "label_origin": label_origin,
        "overall_status": status,
    }


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
        print("Paper06 admissibility classification derivation")
        print("=" * 65)
        rows = []
        for topology in TOPOLOGIES:
            for mode in MODES:
                row = derive_row(topology, mode)
                rows.append(row)
                print(
                    f"  {topology}x{mode}: "
                    f"chi={row['chi_from_lc_curvature']} "
                    f"F_shell={row['F_on_shell']} "
                    f"residual={row['F_on_shell_minus_chi']} "
                    f"label={row['admissibility_label']} "
                    f"status={row['overall_status']}"
                )

        print_rows(rows)

        fail_count = sum(1 for row in rows if row["overall_status"] != "PASS")
        label_counts: dict[str, int] = {}
        for row in rows:
            label = row["admissibility_label"]
            label_counts[label] = label_counts.get(label, 0) + 1

        print()
        print(f"total_branches={len(rows)}")
        print(f"label_counts={label_counts}")
        print(f"fail_count={fail_count}")
        print(f"overall_verdict={'PASS' if fail_count == 0 else 'FAIL'}")
        return 0 if fail_count == 0 else 1
    finally:
        teardown_log()


if __name__ == "__main__":
    raise SystemExit(main())
