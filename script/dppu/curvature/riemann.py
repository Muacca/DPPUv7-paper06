"""
Riemann Tensor Computation and Verification
============================================

Provides Riemann tensor computation in frame basis and strict
antisymmetry verification (3-level approach).

Antisymmetry Properties:
    R_{abcd} must satisfy:
    - R_{abcd} = -R_{bacd}  (antisymmetric in first pair)
    - R_{abcd} = -R_{abdc}  (antisymmetric in second pair)
    - R_{abcd} = R_{cdab}   (pair exchange symmetry)

Verification Strategy (3-level):
    1. Symbolic: Try to prove residual = 0 using simplification
    2. Numerical: Try to find witness point where residual ≠ 0
    3. Classify: PROVED_ZERO, WITNESS_NONZERO, or UNPROVED

Author: Muacca
"""

from typing import Dict, List, Optional, Any

from sympy import S, cancel
from sympy.tensor.array import MutableDenseNDimArray

from ..utils.symbolic import prove_zero, find_nonzero_witness


class RiemannAntisymmetryError(Exception):
    """
    Exception raised when Riemann tensor violates antisymmetry.

    This is a critical error indicating either:
    - Implementation bug in connection/curvature computation
    - Incorrect structure constants for the manifold
    - Numerical instability (rare with symbolic computation)

    Attributes:
        violation_type: Type of violation ('WITNESS_NONZERO' or 'UNPROVED')
        violations: List of violation details

    Example:
        >>> raise RiemannAntisymmetryError(
        ...     "WITNESS_NONZERO",
        ...     [{'indices': (0,1,2,3), 'antisymmetry': 'ab', 'residual': r}]
        ... )
    """

    def __init__(self, violation_type: str, violations: List[Dict]):
        self.violation_type = violation_type
        self.violations = violations

        # Build detailed error message
        msg_lines = ["Riemann tensor antisymmetry violated", ""]
        msg_lines.append(f"Violation Type: {violation_type}")

        if violations:
            msg_lines.append("Violated components (first 3):")
            for v in violations[:3]:
                indices = v['indices']
                residual = v['residual']

                msg_lines.append(
                    f"  [{indices[0]},{indices[1]},{indices[2]},{indices[3]}]: "
                    f"residual = {residual}"
                )

                if 'witness_point' in v:
                    witness = v['witness_point']
                    value = v['witness_value']
                    point_str = ', '.join([f"{k}={val:.3f}" for k, val in witness.items()])
                    msg_lines.append(f"    (witness at {point_str} → value={value:.3e})")

            msg_lines.append("")
            msg_lines.append(f"Total violations: {len(violations)}")

        super().__init__('\n'.join(msg_lines))


def verify_antisymmetry_strict(
    R_abcd: MutableDenseNDimArray,
    dim: int,
    symbols_list: Optional[List] = None,
    logger: Optional[Any] = None
) -> None:
    """
    Strictly verify Riemann antisymmetry with 3-level judgment.

    Checks both antisymmetries:
    - ab: R_{abcd} = -R_{bacd}
    - cd: R_{abcd} = -R_{abdc}

    Args:
        R_abcd: Riemann tensor array (lowered indices)
        dim: Dimension (typically 4)
        symbols_list: (Deprecated) Unused, kept for compatibility
        logger: Optional logger for progress messages

    Raises:
        RiemannAntisymmetryError: If any violation found
    """
    if logger:
        logger.info("Verifying Riemann antisymmetry (STRICT 3-level mode)...")

    violations = []

    # Check ab antisymmetry: R_abcd = -R_bacd
    for a in range(dim):
        for b in range(a + 1, dim):
            for c in range(dim):
                for d in range(dim):
                    residual = R_abcd[a, b, c, d] + R_abcd[b, a, c, d]

                    # Stage 1: Try symbolic proof
                    proved, method = prove_zero(residual)
                    if proved:
                        continue  # PASSED

                    # Stage 2: Try numerical witness
                    residual_symbols = sorted(residual.free_symbols, key=str)
                    has_witness, witness, value = find_nonzero_witness(
                        residual, residual_symbols
                    )

                    if has_witness:
                        violations.append({
                            'type': 'WITNESS_NONZERO',
                            'indices': (a, b, c, d),
                            'antisymmetry': 'ab',
                            'residual': residual,
                            'witness_point': witness,
                            'witness_value': value
                        })
                    else:
                        violations.append({
                            'type': 'UNPROVED',
                            'indices': (a, b, c, d),
                            'antisymmetry': 'ab',
                            'residual': residual
                        })

    # Check cd antisymmetry: R_abcd = -R_abdc
    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(c + 1, dim):
                    residual = R_abcd[a, b, c, d] + R_abcd[a, b, d, c]

                    proved, method = prove_zero(residual)
                    if proved:
                        continue

                    residual_symbols = sorted(residual.free_symbols, key=str)
                    has_witness, witness, value = find_nonzero_witness(
                        residual, residual_symbols
                    )

                    if has_witness:
                        violations.append({
                            'type': 'WITNESS_NONZERO',
                            'indices': (a, b, c, d),
                            'antisymmetry': 'cd',
                            'residual': residual,
                            'witness_point': witness,
                            'witness_value': value
                        })
                    else:
                        violations.append({
                            'type': 'UNPROVED',
                            'indices': (a, b, c, d),
                            'antisymmetry': 'cd',
                            'residual': residual
                        })

    # Report results
    if violations:
        violation_type = violations[0]['type']
        raise RiemannAntisymmetryError(violation_type, violations)
    else:
        if logger:
            logger.success("Riemann antisymmetry: PASSED (all PROVED_ZERO)")


def compute_riemann_tensor(
    Gamma: MutableDenseNDimArray,
    C: MutableDenseNDimArray,
    dim: int,
    frame_deriv: Optional[Any] = None
) -> MutableDenseNDimArray:
    """
    Compute Riemann tensor R^a_{bcd} from connection and structure constants.

    Formula (frame basis):
        R^a_{bcd} = e_c(Γ^a_{bd}) - e_d(Γ^a_{bc}) + Γ^a_{ec} Γ^e_{bd} - Γ^a_{ed} Γ^e_{bc} + Γ^a_{be} C^e_{cd}

    Optimization (Option S):
        Since C^e_{cd} is antisymmetric in (c,d), the identity
            R^a_{bcd} = -R^a_{bdc}
        holds exactly (proof: swap c↔d in formula, use C^e_{dc} = -C^e_{cd}).
        We compute only c < d components and fill the c > d half by antisymmetry,
        reducing cancel() calls from dim^4 to dim^2 * C(dim,2).

    WARNING — Out-of-Memory risk with off-diagonal shear:
        When enable_offdiag_shear=True and all three parameters q₃, q₄, q₅ are
        treated as simultaneous SymPy symbols, each Gamma component expands to
        800–1800 ops as a rational function in (z₃, z₄, z₅, f, r).
        The products Gamma[a,e,c] * Gamma[e,b,d] then produce intermediate
        polynomials of ~3,200,000 ops, causing exponential heap growth.
        This function will exhaust memory (OOM) before completing, even on
        machines with 64 GB RAM.

        Do NOT run the full pipeline with an S3 off-diagonal shear
        configuration that activates more than one off-diagonal qi simultaneously.
        For off-diagonal Hessian computations, use the single-variable Method 2
        approach (z = exp(q/√2) rational substitution, one qi at a time) or a
        numerical finite-difference pipeline instead.

    Args:
        Gamma: Connection coefficients Γ^a_{bc}
        C: Structure constants C^a_{bc}
        dim: Dimension
        frame_deriv: Optional callable f(index, expr) computing the frame derivative e_{index}(expr)

    Returns:
        Riemann tensor R^a_{bcd} as MutableDenseNDimArray
    """
    Riemann = MutableDenseNDimArray.zeros(dim, dim, dim, dim)

    for a in range(dim):
        for b in range(dim):
            # Compute only independent components with c < d; fill c > d by antisymmetry
            for c in range(dim):
                for d in range(c + 1, dim):
                    term = S.Zero

                    if frame_deriv is not None:
                        term += frame_deriv(c, Gamma[a, b, d]) - frame_deriv(d, Gamma[a, b, c])

                    # Γ^a_{ec} Γ^e_{bd} - Γ^a_{ed} Γ^e_{bc}
                    for e in range(dim):
                        term += Gamma[a, e, c] * Gamma[e, b, d]
                        term -= Gamma[a, e, d] * Gamma[e, b, c]

                    # + Γ^a_{be} C^e_{cd}
                    for e in range(dim):
                        term += Gamma[a, b, e] * C[e, c, d]

                    if term != S.Zero:
                        val = cancel(term)
                        Riemann[a, b, c, d] =  val
                        Riemann[a, b, d, c] = -val  # R^a_{bdc} = -R^a_{bcd}

    return Riemann


def lower_first_index(
    R_a_bcd: MutableDenseNDimArray,
    metric: Any,
    dim: int
) -> MutableDenseNDimArray:
    """
    Lower the first index of Riemann tensor: R_{abcd} = η_{ae} R^e_{bcd}

    Args:
        R_a_bcd: Riemann tensor with upper first index
        metric: Frame metric η_{ab}
        dim: Dimension

    Returns:
        Riemann tensor R_{abcd} with all lower indices
    """
    R_abcd = MutableDenseNDimArray.zeros(dim, dim, dim, dim)

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                for d in range(dim):
                    val = S.Zero
                    for e in range(dim):
                        val += metric[a, e] * R_a_bcd[e, b, c, d]
                    R_abcd[a, b, c, d] = val

    return R_abcd
