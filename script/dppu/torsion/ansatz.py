"""
Torsion Ansatz Construction
===========================

Constructs the torsion tensor T_{abc} based on the chosen mode.

Physical Background:
    The torsion tensor can be decomposed into three irreducible parts:
    1. T1 (Axial/Antisymmetric): Totally antisymmetric part
       - Parameterized by axial vector S^μ
       - T_{abc} = (2η/r) ε_{abc} for spatial indices

    2. T2 (Vector-Trace): Trace part
       - Parameterized by vector V_μ
       - T_{abc} = (1/3)(η_{ac} V_b - η_{ab} V_c)

    3. T3 (Tensor): Traceless hook-symmetric part (not used in our ansatz)

Our M³×S¹ Ansatz:
    - T1: S^μ = (η/r)(0, 0, 0, 1) in frame basis
    - T2: V_μ = V · δ^3_μ (only τ-component)

Author: Muacca
"""

from typing import Any, Optional

from sympy import S, Rational, cancel
from sympy.tensor.array import MutableDenseNDimArray

from .mode import Mode
from ..utils.epsilon import epsilon_3d


def construct_torsion_tensor(
    mode: Mode,
    r: Any,
    eta: Any,
    V: Any,
    metric: Any,
    dim: int = 4,
    logger: Optional[Any] = None,
    frame_convention: str = "legacy_euclidean",
    signature: str = "euclidean",
) -> MutableDenseNDimArray:
    """
    Construct torsion tensor T_{abc} based on mode.

    Args:
        mode: Torsion mode (AX, VT, or MX)
        r: Radius parameter (Symbol or number)
        eta: Axial parameter η (used in AX/MX modes)
        V: Vector-trace parameter V (used in VT/MX modes)
        metric: Frame metric η_{ab}
        dim: Dimension (default 4)
        logger: Optional logger
        frame_convention: "legacy_euclidean" or "lz_native"
        signature: "euclidean" or "lorentzian"

    Returns:
        Torsion tensor T_{abc} as MutableDenseNDimArray

    Notes:
        - For AX mode: Only T1 component is non-zero
        - For VT mode: Only T2 component is non-zero
        - For MX mode: Both T1 and T2 are non-zero
        - For signature="lorentzian" + frame_convention="lz_native":
            AX spatial block is {1,2,3}; VT trace vector is at time index 0.
    """
    if signature == "lorentzian" and frame_convention == "lz_native":
        return _construct_torsion_tensor_lz_native(mode, r, eta, V, metric, dim, logger)

    if logger:
        logger.info(f"Constructing Torsion T_abc (Mode: {mode.value})...")

    T_tensor = MutableDenseNDimArray.zeros(dim, dim, dim)

    # T1: Totally Antisymmetric (Axial)
    if mode.has_axial:
        if eta != S.Zero:
            if logger:
                logger.info("  Adding T1 (Axial/Antisymmetric)...")
            T_tensor = _add_axial_torsion(T_tensor, r, eta, dim)

    # T2: Vector Trace
    if mode.has_vector_trace:
        if V != S.Zero:
            if logger:
                logger.info("  Adding T2 (Vector Trace)...")
            T_tensor = _add_vector_trace_torsion(T_tensor, V, metric, dim)

    if logger:
        nonzero = sum(
            1 for a in range(dim)
            for b in range(dim)
            for c in range(dim)
            if T_tensor[a, b, c] != S.Zero
        )
        logger.info(f"  Non-zero torsion components: {nonzero}")
        logger.success("Torsion constructed")

    return T_tensor


def _add_axial_torsion(
    T: MutableDenseNDimArray,
    r: Any,
    eta: Any,
    dim: int
) -> MutableDenseNDimArray:
    """
    Add axial (T1) torsion component (legacy_euclidean convention).

    Formula for M³×S¹:
        T_{ijk} = (2η/r) ε_{ijk} for i,j,k ∈ {0,1,2}

    This corresponds to axial vector S^μ = (0, 0, 0, η) in τ-direction.
    """
    for i in range(3):
        for j in range(3):
            for k in range(3):
                eps = epsilon_3d(i, j, k)
                if eps != 0:
                    T[i, j, k] = T[i, j, k] + 2 * eta * eps / r

    return T


def _add_vector_trace_torsion(
    T: MutableDenseNDimArray,
    V: Any,
    metric: Any,
    dim: int
) -> MutableDenseNDimArray:
    """
    Add vector-trace (T2) torsion component.

    Formula:
        T_{abc} = (1/3)(η_{ac} V_b - η_{ab} V_c)

    Where V_μ = V · δ^3_μ (only τ-component is non-zero).
    """
    # T_vec has only the τ-component non-zero
    T_vec = [S.Zero, S.Zero, S.Zero, V]

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                val = (metric[a, c] * T_vec[b] - metric[a, b] * T_vec[c]) * Rational(1, 3)
                if val != S.Zero:
                    T[a, b, c] = T[a, b, c] + val

    return T


def extract_torsion_parameters(
    T: MutableDenseNDimArray,
    r: Any,
    dim: int = 4
) -> dict:
    """
    Extract torsion parameters (η, V) from a torsion tensor.

    This is the inverse of construct_torsion_tensor.

    Args:
        T: Torsion tensor T_{abc}
        r: Radius parameter
        dim: Dimension

    Returns:
        Dictionary with 'eta' and 'V' parameters
    """
    # Extract η from T_{012} (axial part)
    if T[0, 1, 2] != S.Zero:
        eta = cancel(T[0, 1, 2] * r / 2)
    else:
        eta = S.Zero

    # Extract V from T_{303} (vector-trace part)
    # T_{303} = (1/3)(η_{33} V_0 - η_{30} V_3) = (1/3)(V_0 - 0) for τ-direction
    # Actually, for V_μ = (0,0,0,V), the pattern is different
    # T_{a3c} - T_{a3c} contributions...
    # For simplicity, check T_{030} = (1/3)(η_{00}V_3 - η_{03}V_0) = V/3
    if T[0, 3, 0] != S.Zero:
        V_param = cancel(3 * T[0, 3, 0])
    else:
        V_param = S.Zero

    return {'eta': eta, 'V': V_param}


# ---------------------------------------------------------------------------
# LZ-native torsion ansatz (signature="lorentzian", frame_convention="lz_native")
# ---------------------------------------------------------------------------

def _construct_torsion_tensor_lz_native(
    mode: Mode,
    r: Any,
    eta: Any,
    V: Any,
    metric: Any,
    dim: int,
    logger: Optional[Any] = None,
) -> MutableDenseNDimArray:
    """
    Lorentzian DPPU-LZ native torsion ansatz.

    Index convention:
        time index   = 0
        spatial indices = {1, 2, 3}

    AX:  T_{ijk} = (2η/r) ε_{i-1,j-1,k-1}  for i,j,k ∈ {1,2,3}
    VT:  T_{abc} = (1/3)(g_{ac} V_b - g_{ab} V_c)  with V_a = V δ^0_a
    """
    if logger:
        logger.info(f"Constructing Torsion T_abc [LZ-native] (Mode: {mode.value})...")

    T_tensor = MutableDenseNDimArray.zeros(dim, dim, dim)

    if mode.has_axial and eta != S.Zero:
        if logger:
            logger.info("  Adding T1 (Axial, LZ-native spatial {1,2,3})...")
        T_tensor = _add_axial_torsion_lz_native(T_tensor, r, eta, dim)

    if mode.has_vector_trace and V != S.Zero:
        if logger:
            logger.info("  Adding T2 (Vector Trace, LZ-native time index 0)...")
        T_tensor = _add_vector_trace_torsion_lz_native(T_tensor, V, metric, dim)

    if logger:
        nonzero = sum(
            1 for a in range(dim)
            for b in range(dim)
            for c in range(dim)
            if T_tensor[a, b, c] != S.Zero
        )
        logger.info(f"  Non-zero torsion components: {nonzero}")
        logger.success("Torsion constructed [LZ-native]")

    return T_tensor


def _add_axial_torsion_lz_native(
    T: MutableDenseNDimArray,
    r: Any,
    eta: Any,
    dim: int,
) -> MutableDenseNDimArray:
    """
    Add axial (T1) torsion component in LZ-native convention.

    Formula:
        T_{ijk} = (2η/r) ε_{i-1,j-1,k-1}  for i,j,k ∈ {1,2,3}

    Spatial block {1,2,3} maps to epsilon_3d with 0-indexed offset.
    """
    for i in range(1, 4):
        for j in range(1, 4):
            for k in range(1, 4):
                eps = epsilon_3d(i - 1, j - 1, k - 1)
                if eps != 0:
                    T[i, j, k] = T[i, j, k] + 2 * eta * eps / r
    return T


def _add_vector_trace_torsion_lz_native(
    T: MutableDenseNDimArray,
    V: Any,
    metric: Any,
    dim: int,
) -> MutableDenseNDimArray:
    """
    Add vector-trace (T2) torsion component in LZ-native convention.

    Formula:
        T_{abc} = (1/3)(g_{ac} V_b - g_{ab} V_c)

    V_a = V · δ^0_a  (time component only in LZ-native).
    """
    T_vec = [V, S.Zero, S.Zero, S.Zero]

    for a in range(dim):
        for b in range(dim):
            for c in range(dim):
                val = (metric[a, c] * T_vec[b] - metric[a, b] * T_vec[c]) * Rational(1, 3)
                if val != S.Zero:
                    T[a, b, c] = T[a, b, c] + val

    return T


def torsion_trace(
    T: MutableDenseNDimArray,
    dim: int
) -> list:
    """
    Compute the torsion trace vector.

    Definition:
        T_b = T^a_{ab}

    Args:
        T: Torsion tensor T_{abc}
        dim: Dimension

    Returns:
        Trace vector [T_0, T_1, T_2, T_3]
    """
    trace = []
    for b in range(dim):
        val = S.Zero
        for a in range(dim):
            val += T[a, a, b]
        trace.append(cancel(val))
    return trace
