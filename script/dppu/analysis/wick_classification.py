"""
Wick Classification
===================

Phase 1E: Wick classification schema for DPPU-LZ quantities.

Defines classification labels, confidence levels, and builds the Phase 1E
quantity-level classification table for the 20 mandatory quantities.

Also provides Pontryagin diagnostic case runners for Part D decomposition:
  - Case E0: Legacy Euclidean baseline (reference)
  - Case LZ: Final Lorentzian EC/Hodge result
  - Case H1: LZ curvature + legacy-like combinatorial Hodge (diagnostic only)
  - Case S1: LZ curvature + Euclidean metric-aware Hodge (diagnostic only)
  - Case T1: Legacy torsion ansatz under Lorentzian metric-aware P (diagnostic only)
  - Case I1: Legacy-index remapped Riemann vs native LZ Riemann (diagnostic only)

IMPORTANT:
  - Physical interpretations (P=0 protection, stability, ghost-free) are
    NOT provided in this module. This module classifies algebraic behavior only.
  - Diagnostic cases are marked DIAGNOSTIC_ONLY and must not be treated as
    physical models.
"""

from sympy import Matrix, Rational, S, cancel, factor, symbols
from sympy.tensor.array import MutableDenseNDimArray

from ..utils.epsilon import epsilon_symbol, epsilon_tensor_down


# ---------------------------------------------------------------------------
# Classification labels
# ---------------------------------------------------------------------------

CLASSIFICATION_LABELS = [
    "INVARIANT",
    "SIGN_FLIPPED",
    "I_WEIGHTED",
    "METRIC_DEPENDENT",
    "HODGE_DEPENDENT",
    "INDEX_CONVENTION_DEPENDENT",
    "TORSION_ANSATZ_DEPENDENT",
    "DEFINITION_DEPENDENT",
    "BOUNDARY_SHIFTED",
    "ZERO_TO_ZERO",
    "NONZERO_TO_ZERO",
    "ZERO_TO_NONZERO",
    "NON_CORRESPONDING",
    "UNCLASSIFIED",
]

CONFIDENCE_LABELS = [
    "PROVEN_SYMBOLIC",
    "CHECKED_BY_SCRIPT",
    "DIAGNOSTIC_ONLY",
    "HEURISTIC",
    "UNRESOLVED",
]


# ---------------------------------------------------------------------------
# Pontryagin diagnostic helpers
# ---------------------------------------------------------------------------

def pontryagin_legacy_flat_hodge(R_abcd, metric_inv, dim: int = 4,
                                  simplify_result: bool = True):
    """
    H1 diagnostic: compute P using legacy-style flat index contraction.

    Differs from metric-aware P in that the second curvature factor R_ab^{ef}
    uses NO metric raising on the form indices (ef), i.e., treated as R_{abef}.

    P_flat = (1/8) * sum_{a,b,c,d,e,f} R^{abcd} * eps_symbol(c,d,e,f) * R_{abef}

    where R^{abcd} is still fully raised with metric_inv (Lorentzian).

    For orthonormal frame eps_symbol = eps_tensor_down, so the ONLY difference
    from P_LZ is that the second factor R_ab^{ef} -> R_{abef} (no ef-raise).

    This isolates the contribution of metric raising on the form indices of
    the second curvature factor.

    DIAGNOSTIC ONLY / NOT A PHYSICAL MODEL.
    """
    R_up = _raise_all_indices_local(R_abcd, metric_inv, dim)

    P = S.Zero
    for a in range(dim):
        for b in range(dim):
            if a == b:
                continue
            for c in range(dim):
                for d in range(dim):
                    Rup_abcd = R_up[a, b, c, d]
                    if Rup_abcd == S.Zero:
                        continue
                    for e in range(dim):
                        for f in range(dim):
                            eps = epsilon_symbol(c, d, e, f)
                            if eps == 0:
                                continue
                            R_abef = R_abcd[a, b, e, f]
                            if R_abef != S.Zero:
                                P += Rup_abcd * eps * R_abef

    P = Rational(1, 8) * P
    if simplify_result:
        P = cancel(P)
    return P


def pontryagin_euclidean_metric_aware(R_abcd, dim: int = 4,
                                       simplify_result: bool = True):
    """
    S1 diagnostic: compute metric-aware P but using Euclidean metric (I).

    Uses hodge_dual_2form with the Euclidean identity metric delta_{ab}.
    The curvature R_abcd is taken from the Lorentzian engine.

    This isolates whether the Lorentzian signature (eta_00 = -1) is responsible
    for P_LZ = 0 by comparing:
      P_S1_euclidean = P with delta metric in Hodge (on Lorentzian curvature)
      P_S1_lorentzian = P_LZ = EXACT_ZERO (on same Lorentzian curvature)

    DIAGNOSTIC ONLY / NOT A PHYSICAL MODEL.
    """
    from ..curvature.pontryagin_lz import pontryagin_density_lorentzian

    DELTA = Matrix.eye(dim)
    return pontryagin_density_lorentzian(
        R_abcd, DELTA, simplify_result=simplify_result
    )


def _raise_all_indices_local(R_abcd, metric_inv, dim: int = 4):
    """Raise all four indices of R_abcd using metric_inv."""
    R_up = MutableDenseNDimArray.zeros(dim, dim, dim, dim)
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    val = S.Zero
                    for A in range(dim):
                        gAa = metric_inv[a, A]
                        if gAa == S.Zero:
                            continue
                        for B in range(dim):
                            gBb = metric_inv[b, B]
                            if gBb == S.Zero:
                                continue
                            for C in range(dim):
                                gCc = metric_inv[c, C]
                                if gCc == S.Zero:
                                    continue
                                for D in range(dim):
                                    gDd = metric_inv[d, D]
                                    if gDd == S.Zero:
                                        continue
                                    rv = R_abcd[A, B, C, D]
                                    if rv != S.Zero:
                                        val += gAa * gBb * gCc * gDd * rv
                    if val != S.Zero:
                        R_up[a, b, c, d] = val
    return R_up



# ---------------------------------------------------------------------------
# Quantity-level classification table
# ---------------------------------------------------------------------------

def classify_quantity_wick_behavior(quantity_name: str, evidence: dict) -> dict:
    """
    Return classification labels and confidence for a named quantity.

    Parameters
    ----------
    quantity_name : str
        Name of the quantity (must match a known entry in the table).
    evidence : dict
        Optional additional evidence to merge.

    Returns
    -------
    dict with keys: quantity, old_euclidean, lz_lorentzian, labels,
                    confidence, evidence_ref, notes
    """
    table = {row["quantity"]: row for row in _RAW_CLASSIFICATION_TABLE}
    if quantity_name in table:
        row = dict(table[quantity_name])
        if evidence:
            row["evidence_ref"] = row.get("evidence_ref", "") + " | " + str(evidence)
        return row
    return {
        "quantity": quantity_name,
        "old_euclidean": "unknown",
        "lz_lorentzian": "unknown",
        "labels": ["UNCLASSIFIED"],
        "confidence": "UNRESOLVED",
        "evidence_ref": str(evidence),
        "notes": "quantity not found in classification table",
    }


def build_phase1e_classification_table() -> list:
    """
    Build Wick classification table for Phase 1E quantities.

    Returns list of dicts, one per quantity, with fields:
      quantity, old_euclidean, lz_lorentzian, labels, confidence, evidence_ref, notes
    """
    return list(_RAW_CLASSIFICATION_TABLE)


def compare_pontryagin_diagnostic_cases(cases: list) -> list:
    """
    Run or summarize diagnostic Pontryagin comparison cases.

    Each input case dict specifies: case_id, setup, and optionally
    R_abcd / metric data for live computation.

    Returns list of result dicts with fields:
      case_id, setup, AX, VT, MX, classification, confidence, comment, is_live
    """
    results = []
    for case in cases:
        cid = case.get("case_id", "UNKNOWN")
        if cid in _STATIC_DIAGNOSTIC_CASES:
            results.append(dict(_STATIC_DIAGNOSTIC_CASES[cid]))
        else:
            results.append({
                "case_id": cid,
                "setup": case.get("setup", ""),
                "AX": "unresolved",
                "VT": "unresolved",
                "MX": "unresolved",
                "classification": "UNCLASSIFIED",
                "confidence": "UNRESOLVED",
                "comment": "case not recognized",
                "is_live": False,
            })
    return results


def render_classification_table_markdown(rows: list) -> str:
    """
    Render Wick classification table as Markdown.

    Produces a GitHub-Flavored Markdown table with columns:
    Quantity | Old Euclidean | DPPU-LZ Lorentzian | Labels | Confidence | Evidence | Notes
    """
    header = (
        "| Quantity | Old Euclidean behavior | DPPU-LZ Lorentzian behavior "
        "| Labels | Confidence | Evidence | Notes |\n"
        "|---|---|---|---|---|---|---|\n"
    )
    lines = [header]
    for row in rows:
        labels = ", ".join(row.get("labels", []))
        line = (
            f"| {row.get('quantity','?')} "
            f"| {row.get('old_euclidean','?')} "
            f"| {row.get('lz_lorentzian','?')} "
            f"| {labels} "
            f"| {row.get('confidence','?')} "
            f"| {row.get('evidence_ref','?')} "
            f"| {row.get('notes','')} |\n"
        )
        lines.append(line)
    return "".join(lines)


def render_diagnostic_table_markdown(rows: list) -> str:
    """
    Render Pontryagin diagnostic cases table as Markdown.
    """
    header = (
        "| Case ID | Setup | AX | VT | MX "
        "| Classification | Confidence | Comment |\n"
        "|---|---|---|---|---|---|---|---|\n"
    )
    lines = [header]
    for row in rows:
        line = (
            f"| {row.get('case_id','?')} "
            f"| {row.get('setup','?')} "
            f"| {row.get('AX','?')} "
            f"| {row.get('VT','?')} "
            f"| {row.get('MX','?')} "
            f"| {row.get('classification','?')} "
            f"| {row.get('confidence','?')} "
            f"| {row.get('comment','')} |\n"
        )
        lines.append(line)
    return "".join(lines)


# ---------------------------------------------------------------------------
# Raw classification data (20 mandatory quantities)
# ---------------------------------------------------------------------------

_RAW_CLASSIFICATION_TABLE = [
    {
        "quantity": "frame metric eta_ab",
        "old_euclidean": "delta_ab = diag(+1,+1,+1,+1)",
        "lz_lorentzian": "eta_ab = diag(-1,+1,+1,+1)",
        "labels": ["SIGN_FLIPPED", "METRIC_DEPENDENT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1B, Phase1C",
        "notes": "eta_00 flipped from +1 to -1; root cause of Lorentzian sign differences",
    },
    {
        "quantity": "coframe e^a",
        "old_euclidean": "e^0,e^1,e^2 spatial; e^3 fiber/time (old index 3=fiber)",
        "lz_lorentzian": "e^0 time; e^1,e^2,e^3 spatial (LZ index 0=time)",
        "labels": ["INDEX_CONVENTION_DEPENDENT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1A",
        "notes": "index relabeling: a_new = (a_old+1) mod 4; orientation sign recorded in Phase1A",
    },
    {
        "quantity": "volume/lapse factor L vs N(t)",
        "old_euclidean": "lapse L (Euclidean)",
        "lz_lorentzian": "lapse N(t) (Lorentzian, carries -1 from g_00=-N^2)",
        "labels": ["METRIC_DEPENDENT", "DEFINITION_DEPENDENT"],
        "confidence": "HEURISTIC",
        "evidence_ref": "Phase1C smoke test (EH/FLRW)",
        "notes": "Lorentzian lapse sign change affects action density; Phase 2 target for full derivation",
    },
    {
        "quantity": "epsilon symbol eps~_abcd",
        "old_euclidean": "eps~_0123 = +1 (combinatorial, signature-independent)",
        "lz_lorentzian": "eps~_0123 = +1 (unchanged)",
        "labels": ["INVARIANT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1B",
        "notes": "pure combinatorial Levi-Civita symbol; no metric dependence",
    },
    {
        "quantity": "epsilon tensor down eps_abcd",
        "old_euclidean": "eps_0123 = +1 (orthonormal frame, sqrt|det g|=1)",
        "lz_lorentzian": "eps_0123 = +1 (orthonormal frame, |det eta|=1)",
        "labels": ["INVARIANT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1B (epsilon_tensor_down(0,1,2,3)=+1 confirmed)",
        "notes": "for orthonormal frame |det g|=1, eps_tensor_down = eps_symbol",
    },
    {
        "quantity": "epsilon tensor up eps^abcd",
        "old_euclidean": "eps^0123 = +1 (delta^{ab} raising: all +1)",
        "lz_lorentzian": "eps^0123 = -1 (eta^{ab} raising: eta^{00}=-1)",
        "labels": ["SIGN_FLIPPED", "METRIC_DEPENDENT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1B (epsilon_tensor_up(0,1,2,3)=-1 confirmed)",
        "notes": "sign flip due to eta^{00}=-1; independent of legacy permutation sign",
    },
    {
        "quantity": "Hodge dual on 2-forms",
        "old_euclidean": "compute_hodge_dual: combinatorial eps_symbol, no metric raising",
        "lz_lorentzian": "hodge_dual_2form: metric-aware, uses eps_tensor_down",
        "labels": ["METRIC_DEPENDENT", "HODGE_DEPENDENT", "DEFINITION_DEPENDENT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1B (Hodge basis table confirmed)",
        "notes": "definition change from combinatorial to metric-aware is key change in Pontryagin",
    },
    {
        "quantity": "double Hodge ** on 2-forms",
        "old_euclidean": "**F = +F  (**=+1, Euclidean 4D)",
        "lz_lorentzian": "**F = -F  (**=-1, Lorentzian 4D)",
        "labels": ["SIGN_FLIPPED", "METRIC_DEPENDENT", "HODGE_DEPENDENT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1B (all 6 basis 2-forms: **F=-F confirmed)",
        "notes": "direct consequence of Lorentzian signature; affects SD/ASD decomposition",
    },
    {
        "quantity": "structure constants C^a_bc",
        "old_euclidean": "S3: spatial block {0,1,2}; Nil3: (2,0,1),(2,1,0)",
        "lz_lorentzian": "S3: spatial block {1,2,3}; Nil3: (3,1,2),(3,2,1)",
        "labels": ["INDEX_CONVENTION_DEPENDENT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1C-Repair (LZ-native remap confirmed for S3, Nil3)",
        "notes": "T3 unaffected (C=0); remap via a_new=(a_old+1) mod 4",
    },
    {
        "quantity": "torsion ansatz T_abc",
        "old_euclidean": "AX: spatial {0,1,2}; VT: vector at index 3 (fiber)",
        "lz_lorentzian": "AX: spatial {1,2,3}; VT: vector at index 0 (time)",
        "labels": ["TORSION_ANSATZ_DEPENDENT", "INDEX_CONVENTION_DEPENDENT",
                   "METRIC_DEPENDENT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1C-Repair (K+K=0 antisymmetry restored after LZ-native repair)",
        "notes": (
            "Lorentzian VT has V at time index 0; "
            "metric lowering sign change for time component was root cause of B2 blocker"
        ),
    },
    {
        "quantity": "contortion K^a_bc / K_abc",
        "old_euclidean": "K_abc without metric raising (identity raise, Euclidean)",
        "lz_lorentzian": "K^a_bc = g^{ad} K_dbc (metric-aware raising with eta)",
        "labels": ["METRIC_DEPENDENT", "TORSION_ANSATZ_DEPENDENT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1C-Repair (K_abc+K_bac=0 confirmed for AX/VT/MX)",
        "notes": "metric-aware raising required for so(1,3) antisymmetry condition",
    },
    {
        "quantity": "EC connection Gamma_EC",
        "old_euclidean": "Gamma_EC = Gamma_LC + K  (Euclidean metric, legacy torsion)",
        "lz_lorentzian": "Gamma_EC = Gamma_LC + K  (Lorentzian metric, LZ-native torsion)",
        "labels": ["METRIC_DEPENDENT", "TORSION_ANSATZ_DEPENDENT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1C-Repair (Gamma_EC,abc+Gamma_EC,bac=0 confirmed)",
        "notes": "so(1,3) metricity condition Gamma_EC,abc+Gamma_EC,bac=0 restored after repair",
    },
    {
        "quantity": "Riemann tensor R_abcd",
        "old_euclidean": "R_abcd (Euclidean metric, legacy structure constants, Euclidean K)",
        "lz_lorentzian": "R_abcd (Lorentzian metric, LZ-native C, Lorentzian K)",
        "labels": ["METRIC_DEPENDENT", "TORSION_ANSATZ_DEPENDENT",
                   "INDEX_CONVENTION_DEPENDENT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1C-Repair (R_abcd+R_bacd=0 confirmed for all modes)",
        "notes": "antisymmetry R_abcd+R_bacd=0 confirmed; R_abcd+R_abdc=0 confirmed",
    },
    {
        "quantity": "Ricci scalar R",
        "old_euclidean": "R (Euclidean)",
        "lz_lorentzian": "R (Lorentzian)",
        "labels": ["METRIC_DEPENDENT", "TORSION_ANSATZ_DEPENDENT"],
        "confidence": "HEURISTIC",
        "evidence_ref": "Phase1C smoke test (flat T3: R=0 confirmed)",
        "notes": "full symbolic evaluation deferred to Phase 2; flat limit verified",
    },
    {
        "quantity": "Weyl scalar C^2",
        "old_euclidean": "C^2 (Euclidean)",
        "lz_lorentzian": "C^2 (Lorentzian)",
        "labels": ["METRIC_DEPENDENT"],
        "confidence": "HEURISTIC",
        "evidence_ref": "Phase1C smoke test (flat T3: C^2=0 confirmed)",
        "notes": "full evaluation deferred; flat limit verified",
    },
    {
        "quantity": "Nieh-Yan density N",
        "old_euclidean": "N (Euclidean, legacy index/torsion)",
        "lz_lorentzian": "N (Lorentzian, LZ-native)",
        "labels": ["INDEX_CONVENTION_DEPENDENT", "TORSION_ANSATZ_DEPENDENT",
                   "METRIC_DEPENDENT"],
        "confidence": "HEURISTIC",
        "evidence_ref": "Phase1C smoke test (flat T3: N=0 confirmed)",
        "notes": "boundary term structure; Lorentzian behavior partially unresolved; Phase 2 target",
    },
    {
        "quantity": "Pontryagin density P_EC",
        "old_euclidean": "AX=zero, VT=zero, MX=nonzero  (Euclidean DPPU report)",
        "lz_lorentzian": (
            "AX=EXACT_ZERO, VT=EXACT_ZERO, MX=EXACT_ZERO  "
            "(under adopted Lorentzian EC/Hodge definition)"
        ),
        "labels": [
            "NONZERO_TO_ZERO",      # for MX
            "ZERO_TO_ZERO",         # for AX, VT
            "HODGE_DEPENDENT",
            "METRIC_DEPENDENT",
            "TORSION_ANSATZ_DEPENDENT",
            "DEFINITION_DEPENDENT",
        ],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": (
            "Phase1D (9/9 EXACT_ZERO confirmed); "
            "Phase1D report S10 (old Euclidean comparison)"
        ),
        "notes": (
            "MX: NONZERO_TO_ZERO under adopted Lorentzian EC/Hodge definition. "
            "Likely source: combination of metric-aware Hodge, Lorentzian signature, "
            "and LZ-native so(1,3) contortion repair. "
            "Further isolation requires diagnostic cases H1/S1/T1. "
            "Physical interpretation NOT provided (Phase 2 target)."
        ),
    },
    {
        "quantity": "self-duality relation",
        "old_euclidean": "omega+= (1/2)(omega + *omega), real eigenvalues +/-1 (**=+1)",
        "lz_lorentzian": "omega+= (1/2)(omega - i*omega), complex eigenvalues +/-i (**=-1)",
        "labels": ["SIGN_FLIPPED", "I_WEIGHTED", "HODGE_DEPENDENT"],
        "confidence": "PROVEN_SYMBOLIC",
        "evidence_ref": "Phase1B (A2-lite definition); conventions v1.1",
        "notes": "Lorentzian SD/ASD is complex-valued; ** eigenvalues are +/-i not +/-1",
    },
    {
        "quantity": "action / reduced potential",
        "old_euclidean": "S_E (Euclidean action, positive definite signature)",
        "lz_lorentzian": "S_L (Lorentzian action, lapse N carries -1)",
        "labels": ["METRIC_DEPENDENT", "SIGN_FLIPPED"],
        "confidence": "HEURISTIC",
        "evidence_ref": "Phase1C smoke test (EH/FLRW sign structure verified)",
        "notes": "V_eff sign reinterpreted under Lorentzian FLRW ansatz; Phase 2 target",
    },
    {
        "quantity": "Hamiltonian constraint status",
        "old_euclidean": "not evaluated here",
        "lz_lorentzian": "not evaluated / Phase 2 target",
        "labels": ["NON_CORRESPONDING"],
        "confidence": "UNRESOLVED",
        "evidence_ref": "Phase1 plan (explicitly deferred to Phase 2)",
        "notes": "Hamiltonian constraint derivation is Phase 2 target; NOT evaluated in Phase 1",
    },
]


# ---------------------------------------------------------------------------
# Static diagnostic case records
# ---------------------------------------------------------------------------

_STATIC_DIAGNOSTIC_CASES = {
    "E0": {
        "case_id": "E0",
        "setup": (
            "Legacy Euclidean baseline: "
            "signature=euclidean, frame_convention=legacy_euclidean, "
            "torsion=legacy, Hodge=compute_hodge_dual (combinatorial epsilon_symbol)"
        ),
        "AX": "zero",
        "VT": "zero",
        "MX": "nonzero",
        "classification": (
            "old reference; AX/VT: zero (consistent with LZ result); "
            "MX: nonzero (changes to EXACT_ZERO in LZ)"
        ),
        "confidence": "CHECKED_BY_SCRIPT",
        "comment": (
            "Result from old DPPU Euclidean implementation. "
            "Referenced from Phase1D report S10. "
            "Re-verification possible via compute_P_from_riemann + Euclidean engine."
        ),
        "is_live": False,
    },
    "LZ": {
        "case_id": "LZ",
        "setup": (
            "Final Lorentzian result: "
            "signature=lorentzian, frame_convention=lz_native, "
            "torsion=LZ-native so(1,3), Hodge=hodge_dual_2form (metric-aware)"
        ),
        "AX": "EXACT_ZERO",
        "VT": "EXACT_ZERO",
        "MX": "EXACT_ZERO",
        "classification": (
            "AX/VT: ZERO_TO_ZERO (same as E0); "
            "MX: NONZERO_TO_ZERO (differs from E0)"
        ),
        "confidence": "PROVEN_SYMBOLIC",
        "comment": (
            "From Phase1D: 9/9 mandatory cases classified as EXACT_ZERO "
            "under adopted Lorentzian EC/Hodge definition. "
            "cancel() sufficient; factor()/simplify() not needed."
        ),
        "is_live": False,
    },
    "H1": {
        "case_id": "H1",
        "setup": (
            "H1 diagnostic (diagnostic only / not a physical model): "
            "R = LZ Lorentzian EC curvature (T3 MX), "
            "P_metric_hodge = hodge_dual_2form (metric-aware, EXACT_ZERO), "
            "P_legacy_flat = legacy-like flat contraction (no metric raise on ef indices)"
        ),
        "AX": "EXACT_ZERO (metric-aware) | computed live (flat)",
        "VT": "EXACT_ZERO (metric-aware) | computed live (flat)",
        "MX": "EXACT_ZERO (metric-aware) | computed live (flat)",
        "classification": "HODGE_DEPENDENT (pending live computation)",
        "confidence": "DIAGNOSTIC_ONLY",
        "comment": (
            "Live computation performed in check script. "
            "If P_flat != 0 and P_metric_aware = 0: "
            "metric raising on ef indices is the discriminating factor. "
            "diagnostic only / not a physical model."
        ),
        "is_live": True,
    },
    "S1": {
        "case_id": "S1",
        "setup": (
            "S1 diagnostic (diagnostic only / not a physical model): "
            "Same LZ curvature R_abcd (T3 MX); "
            "P_euclidean_metric = hodge_dual_2form with delta metric (I); "
            "P_lorentzian_metric = hodge_dual_2form with eta metric (EXACT_ZERO)"
        ),
        "AX": "EXACT_ZERO (lorentzian) | computed live (euclidean metric)",
        "VT": "EXACT_ZERO (lorentzian) | computed live (euclidean metric)",
        "MX": "EXACT_ZERO (lorentzian) | computed live (euclidean metric)",
        "classification": "METRIC_DEPENDENT / SIGN_FLIPPED (pending live computation)",
        "confidence": "DIAGNOSTIC_ONLY",
        "comment": (
            "If P_S1_euc != 0 and P_S1_lor = 0: "
            "eta_00 = -1 sign in metric-aware Hodge is discriminating factor. "
            "diagnostic only / not a physical model."
        ),
        "is_live": True,
    },
    "T1": {
        "case_id": "T1",
        "setup": (
            "T1 diagnostic (diagnostic only / not a physical model): "
            "Legacy torsion ansatz (spatial {0,1,2}, VT at index 3) "
            "vs LZ-native torsion (spatial {1,2,3}, VT at index 0); "
            "same Lorentzian metric-aware P definition"
        ),
        "AX": "EXACT_ZERO (LZ-native) | computed live (legacy torsion)",
        "VT": "EXACT_ZERO (LZ-native) | computed live (legacy torsion)",
        "MX": "EXACT_ZERO (LZ-native) | computed live (legacy torsion)",
        "classification": "TORSION_ANSATZ_DEPENDENT (pending live computation)",
        "confidence": "DIAGNOSTIC_ONLY",
        "comment": (
            "Legacy torsion violates Lorentzian so(1,3) K+K=0 condition "
            "(confirmed in Phase1C-Audit). "
            "Result with legacy torsion is diagnostic only. "
            "diagnostic only / not a physical model."
        ),
        "is_live": True,
    },
    "I1": {
        "case_id": "I1",
        "setup": (
            "I1 diagnostic (diagnostic only / not a physical model): "
            "Legacy-index Riemann (old {0,1,2,3}) remapped to LZ "
            "vs native LZ Riemann; T3 (trivial C=0), S3 (nontrivial)"
        ),
        "AX": "EXACT_ZERO (native) | computed live (remapped)",
        "VT": "EXACT_ZERO (native) | computed live (remapped)",
        "MX": "EXACT_ZERO (native) | computed live (remapped)",
        "classification": "INDEX_CONVENTION_DEPENDENT (pending live computation)",
        "confidence": "DIAGNOSTIC_ONLY",
        "comment": (
            "T3: trivial (C=0, remap has no effect). "
            "S3: remap changes structure constant keys {0,1,2}->{1,2,3}. "
            "If P_remapped = P_native: index remap alone does not explain P difference. "
            "diagnostic only / not a physical model."
        ),
        "is_live": True,
    },
}
