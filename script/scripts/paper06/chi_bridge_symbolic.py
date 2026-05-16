"""Derive the paper06 chi-bridge identity.

It derives chi from the Levi-Civita spatial curvature of the Appendix A
structure constants, derives the raw MX internal-pair diagnostic from the
bundled EC+NY engine, extracts C from that raw diagnostic without using chi,
and checks C + 9*chi = 0.

Outputs:
  data/chi_bridge_symbolic_YYYYMMDD_HHMMSS.log
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp
from sympy import Matrix, Rational, S
from sympy.tensor.array import MutableDenseNDimArray

HERE = Path(__file__).resolve().parent
SCRIPT_ROOT = HERE.parents[1]     # .../script/
PROJECT_ROOT = HERE.parents[2]    # .../20_DPPUv7-paper06/
DATA_DIR = PROJECT_ROOT / "data"

if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from dppu.curvature.ricci import compute_ricci_scalar_from_tensor, compute_ricci_tensor
from dppu.curvature.riemann import compute_riemann_tensor, lower_first_index
from dppu.engine.contortion import compute_contortion
from dppu.engine.ec_connection import compute_ec_connection
from dppu.engine.levi_civita import compute_christoffel_frame
from dppu.torsion.ansatz import construct_torsion_tensor
from dppu.torsion.mode import Mode
from dppu.topology.lz_invariants import A, B, C, Q, U, W, topology_lz_parameter_table
from dppu.action.reduced_sector import ETA, V_TORSION
from dppu.utils.epsilon import epsilon_tensor_up
from dppu.utils.tee_logger import setup_log, teardown_log


DIM = 4
METRIC = Matrix.diag(-1, 1, 1, 1)
METRIC_INV = METRIC.inv()
TOPOLOGY_PARAMS = topology_lz_parameter_table()

ROW_FIELDS = [
    "row_type", "topology", "a", "b", "c", "u", "v",
    "chi_from_lc_curvature", "P_int_MX_raw", "C_extracted_from_P_int",
    "template_residual", "C_plus_9chi_residual", "overall_status",
]


def build_family_structure_constants() -> MutableDenseNDimArray:
    """Build C^a_{bc} for the five-parameter LZ-native family."""
    arr = MutableDenseNDimArray.zeros(DIM, DIM, DIM)
    for upper, lower1, lower2, value in [
        (1, 2, 3, A),
        (2, 3, 1, B),
        (3, 1, 2, C),
        (2, 1, 2, U),
        (3, 1, 3, W),
    ]:
        arr[upper, lower1, lower2] = value / Q
        arr[upper, lower2, lower1] = -value / Q
    return arr


def compute_internal_pair_pint(R_abcd: MutableDenseNDimArray) -> sp.Expr:
    """Compute the internal-pair P_int diagnostic used in paper06."""
    total = S.Zero
    for c_idx in range(DIM):
        for d_idx in range(c_idx + 1, DIM):
            block_sum = S.Zero
            for a_idx in range(DIM):
                for b_idx in range(DIM):
                    left = R_abcd[a_idx, b_idx, c_idx, d_idx]
                    if left == S.Zero:
                        continue
                    for e_idx in range(DIM):
                        for f_idx in range(DIM):
                            eps = epsilon_tensor_up(METRIC_INV, a_idx, b_idx, e_idx, f_idx)
                            if eps == S.Zero:
                                continue
                            right = R_abcd[e_idx, f_idx, c_idx, d_idx]
                            if right != S.Zero:
                                block_sum += eps * left * right
            if block_sum != S.Zero:
                total += Rational(1, 8) * block_sum
    return sp.factor(sp.cancel(total))


def derive_family_identity() -> dict[str, sp.Expr | str]:
    """Derive chi, P_int, C, and the bridge residual for the full family."""
    structure = build_family_structure_constants()
    gamma_lc = compute_christoffel_frame(structure, DIM, metric=METRIC, metric_inv=METRIC_INV)
    riemann_lc = compute_riemann_tensor(gamma_lc, structure, DIM)
    ricci_lc = compute_ricci_tensor(riemann_lc, DIM)
    scalar_lc = compute_ricci_scalar_from_tensor(ricci_lc, METRIC, DIM)
    chi_expr = sp.factor(sp.cancel(sp.simplify(scalar_lc * Q**2 / 6)))

    torsion = construct_torsion_tensor(
        Mode.MX,
        Q,
        ETA,
        V_TORSION,
        METRIC,
        DIM,
        frame_convention="lz_native",
        signature="lorentzian",
    )
    contortion = compute_contortion(torsion, DIM, metric=METRIC)
    gamma_ec = compute_ec_connection(gamma_lc, contortion, DIM)
    riemann_ec = compute_riemann_tensor(gamma_ec, structure, DIM)
    r_abcd_ec = lower_first_index(riemann_ec, METRIC, DIM)
    pint_expr = compute_internal_pair_pint(r_abcd_ec)

    pint_poly = sp.Poly(sp.expand(pint_expr), V_TORSION, ETA)
    mixed_linear_coeff = pint_poly.coeff_monomial(V_TORSION * ETA)
    c_expr = sp.factor(sp.cancel(sp.simplify(mixed_linear_coeff * (9 * Q**3) / 2)))
    template = sp.factor(
        sp.cancel(
            2 * V_TORSION * ETA
            * (-V_TORSION**2 * Q**2 + 9 * ETA**2 + c_expr)
            / (9 * Q**3)
        )
    )
    template_residual = sp.factor(sp.cancel(sp.simplify(pint_expr - template)))
    bridge_residual = sp.factor(sp.cancel(sp.simplify(c_expr + 9 * chi_expr)))
    status = "PASS" if template_residual == 0 and bridge_residual == 0 else "FAIL"

    return {
        "chi": chi_expr,
        "pint": pint_expr,
        "c_expr": c_expr,
        "template_residual": template_residual,
        "bridge_residual": bridge_residual,
        "status": status,
    }


def make_rows(identity: dict[str, sp.Expr | str]) -> list[dict[str, str]]:
    rows = [
        {
            "row_type": "family",
            "topology": "FIVE_PARAMETER_FAMILY",
            "a": "a",
            "b": "b",
            "c": "c",
            "u": "u",
            "v": "v",
            "chi_from_lc_curvature": str(identity["chi"]),
            "P_int_MX_raw": str(identity["pint"]),
            "C_extracted_from_P_int": str(identity["c_expr"]),
            "template_residual": str(identity["template_residual"]),
            "C_plus_9chi_residual": str(identity["bridge_residual"]),
            "overall_status": str(identity["status"]),
        }
    ]
    chi_expr = identity["chi"]
    c_expr = identity["c_expr"]
    for topology, params in TOPOLOGY_PARAMS.items():
        chi_value = sp.factor(sp.cancel(sp.simplify(chi_expr.subs(params))))
        c_value = sp.factor(sp.cancel(sp.simplify(c_expr.subs(params))))
        residual = sp.factor(sp.cancel(sp.simplify(c_value + 9 * chi_value)))
        status = "PASS" if residual == 0 else "FAIL"
        rows.append(
            {
                "row_type": "specialization",
                "topology": topology,
                "a": str(params[A]),
                "b": str(params[B]),
                "c": str(params[C]),
                "u": str(params[U]),
                "v": str(params[W]),
                "chi_from_lc_curvature": str(chi_value),
                "P_int_MX_raw": "",
                "C_extracted_from_P_int": str(c_value),
                "template_residual": str(identity["template_residual"]),
                "C_plus_9chi_residual": str(residual),
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
        print("Paper06 chi-bridge derivation")
        print("=" * 65)
        identity = derive_family_identity()
        rows = make_rows(identity)
        print_rows(rows)

        fail_count = sum(1 for row in rows if row["overall_status"] != "PASS")
        print(f"family_status={identity['status']}")
        print(f"template_residual={identity['template_residual']}")
        print(f"C_plus_9chi_residual={identity['bridge_residual']}")
        print(f"specialization_count={len(rows) - 1}")
        print(f"fail_count={fail_count}")
        print(f"overall_verdict={'PASS' if fail_count == 0 else 'FAIL'}")
        return 0 if fail_count == 0 else 1
    finally:
        teardown_log()


if __name__ == "__main__":
    raise SystemExit(main())
