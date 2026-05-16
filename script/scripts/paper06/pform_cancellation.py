"""
Paper06 check (Table B1 check 3):
  P_form cancellation: for AX/VT/MX torsion branches across all four topologies
  (S3/T3/Nil3/Sol3), the Pontryagin density P_form = <R,*R> reduces to exact zero.
  This confirms block-diagonal frame orthogonality of the curvature 2-forms.

Engine is used only through step E4.7 (EC Riemann tensor).
Steps E4.3a and E4.3b (LC Weyl) are skipped as they are not needed for P_form.

Outputs:
  data/pform_cancellation_YYYYMMDD_HHMMSS.log  (via tee_logger)
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT_ROOT = HERE.parents[1]     # .../script/
PROJECT_ROOT = HERE.parents[2]    # .../20_DPPUv7-paper06/
DATA_DIR = PROJECT_ROOT / "data"

if str(SCRIPT_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_ROOT))

from dppu.topology.unified import make_engine, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant
from dppu.curvature.pontryagin_lz import (
    pontryagin_density_lorentzian,
    classify_pontryagin_expression,
)
from dppu.utils.tee_logger import setup_log, teardown_log

# Topologies and torsion branches to check
TOPOLOGIES = [
    TopologyType.S3,
    TopologyType.T3,
    TopologyType.NIL3,
    TopologyType.SOL3,
]
TOPO_LABELS = {
    TopologyType.S3:   "S3",
    TopologyType.T3:   "T3",
    TopologyType.NIL3: "Nil3",
    TopologyType.SOL3: "Sol3",
}
TORSION_MODES = [Mode.AX, Mode.VT, Mode.MX]


def run_pipeline_to_riemann(topo: TopologyType, mode: Mode):
    """
    Build engine for the given topology/mode and run only steps E4.1-E4.7.
    Returns (engine, riemann_abcd, metric) after the EC Riemann step.
    Steps E4.3a, E4.3b (LC Weyl) are skipped -- not needed for P_form.
    """
    engine = make_engine(
        topo,
        torsion_mode=mode,
        ny_variant=NyVariant.FULL,
        signature="lorentzian",
        frame_convention="lz_native",
        fiber_mode=FiberMode.NONE,
        enable_squash=False,
    )
    # Run minimum steps to reach EC Riemann tensor
    engine.step_E4_1_setup()
    engine.step_E4_2_metric_and_frame()
    engine.step_E4_3_christoffel_frame()
    # skip E4.3a (LC Riemann/Ricci) and E4.3b (Weyl LC)
    engine.step_E4_4_torsion_ansatz_frame()
    engine.step_E4_5_contortion_frame()
    engine.step_E4_6_ec_connection_frame()
    engine.step_E4_7_riemann_tensor_frame()

    R_abcd = engine.data["riemann_abcd"]
    metric = engine.data["metric_frame"]
    return engine, R_abcd, metric


def main() -> int:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    log_path = setup_log(__file__, log_dir=str(DATA_DIR))
    try:
        print("Paper06 P_form cancellation check (Table B1 check 3)")
        print("=" * 65)

        rows = []
        all_pass = True

        for topo in TOPOLOGIES:
            topo_label = TOPO_LABELS[topo]
            for mode in TORSION_MODES:
                mode_label = mode.value
                tag = f"{topo_label}x{mode_label}"
                print(f"computing: topology={topo_label}  mode={mode_label}")

                try:
                    _, R_abcd, metric = run_pipeline_to_riemann(topo, mode)
                except Exception as exc:
                    print(f"  ERROR during pipeline: {exc}")
                    pform_status = "ENGINE_ERROR"
                    pform_value  = str(exc)
                    row_ok = False
                else:
                    try:
                        P = pontryagin_density_lorentzian(
                            R_abcd, metric, simplify_result=True
                        )
                        pform_status, pform_expr = classify_pontryagin_expression(P)
                        pform_value = str(pform_expr)
                        row_ok = (pform_status == "EXACT_ZERO")
                    except Exception as exc:
                        print(f"  ERROR during P_form computation: {exc}")
                        pform_status = "COMPUTE_ERROR"
                        pform_value  = str(exc)
                        row_ok = False

                if not row_ok:
                    all_pass = False

                overall = "PASS" if row_ok else "FAIL"
                print(
                    f"  pform_status={pform_status}"
                    f"  overall={overall}"
                )
                rows.append({
                    "topology":      topo_label,
                    "torsion_mode":  mode_label,
                    "pform_status":  pform_status,
                    "pform_value":   pform_value,
                    "overall_status": overall,
                })
                print("  row_detail " + json.dumps(rows[-1], sort_keys=True, ensure_ascii=True))

        fail_count = sum(1 for r in rows if r["overall_status"] != "PASS")
        verdict    = "PASS" if all_pass else "FAIL"

        print()
        print(f"total_combinations={len(rows)}")
        print(f"fail_count={fail_count}")
        print(f"overall_verdict={verdict}")

        return 0 if verdict == "PASS" else 1
    finally:
        teardown_log()


if __name__ == "__main__":
    raise SystemExit(main())
