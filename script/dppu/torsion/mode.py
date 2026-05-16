"""
Torsion Mode Enum
=================

Defines the torsion ansatz modes for Einstein-Cartan gravity.

Physical Background:
    The torsion tensor T^λ_μν can be decomposed into three irreducible parts:
    1. Totally antisymmetric (axial/pseudoscalar): S^μ
    2. Vector trace: T_μ
    3. Tensor part (traceless, hook-symmetric)

    For our M³×S¹ ansatz, we consider only the first two:
    - T1 (Axial): S^μ = (η/r)(0,0,0,1) in frame basis
    - T2 (Vector-trace): T_μ = V·δ^3_μ

Modes:
    - AX: Axial-only (S^μ ≠ 0, T_μ = 0)
    - VT: Vector-trace-only (T_μ ≠ 0, S^μ = 0)
    - MX: Mixed (both non-zero)

Author: Muacca
"""

from enum import Enum


class Mode(Enum):
    """
    Torsion ansatz mode enumeration.

    Attributes:
        AX: Axial-only mode.
            Only the totally antisymmetric (axial) component is non-zero.
            Parameters: η ≠ 0, V = 0

        VT: Vector-trace-only mode.
            Only the vector-trace component is non-zero.
            Parameters: η = 0, V ≠ 0

        MX: Mixed mode.
            Both axial and vector-trace components are non-zero.
            Parameters: η ≠ 0, V ≠ 0

    Examples:
        >>> mode = Mode.MX
        >>> mode.value
        'MX'
        >>> Mode['AX']
        <Mode.AX: 'AX'>
    """
    AX = "AX"    # Axial-only (S^μ ≠ 0, T_μ = 0)
    VT = "VT"    # Vector-trace-only (T_μ ≠ 0, S^μ = 0)
    MX = "MX"    # Mixed (both non-zero)

    def __str__(self) -> str:
        return self.value

    @property
    def has_axial(self) -> bool:
        """True if this mode includes axial (T1) component."""
        return self in (Mode.AX, Mode.MX)

    @property
    def has_vector_trace(self) -> bool:
        """True if this mode includes vector-trace (T2) component."""
        return self in (Mode.VT, Mode.MX)

    @property
    def description(self) -> str:
        """Human-readable description of the mode."""
        descriptions = {
            Mode.AX: "Axial-only (T1)",
            Mode.VT: "Vector-trace-only (T2)",
            Mode.MX: "Mixed (T1 + T2)",
        }
        return descriptions[self]


def normalize_mode(label: str) -> str:
    """Normalize a torsion mode label to "AX", "VT", "MX", or "NA".

    Args:
        label: Case-insensitive mode string ("ax", "VT", "MX", "NA", ...).

    Returns:
        Canonical uppercase label.

    Raises:
        ValueError: If *label* does not match any known mode.
    """
    key = label.strip().upper()
    if key == "NA":
        return "NA"
    try:
        return Mode[key].value
    except KeyError:
        raise ValueError(f"unknown torsion mode label: {label!r}")
