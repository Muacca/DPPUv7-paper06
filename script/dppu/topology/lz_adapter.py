"""
LZ-Native Topology Adapter
==========================

Converts old-DPPU structure-constant data (Euclidean frame convention) to
the DPPU-LZ native frame convention (Lorentzian, time index = 0).

Index permutation (Phase 1A result):
    old_to_lz: a_new = (a_old + 1) % 4
    lz_to_old: a_old = (a_new + 3) % 4

Map table:
    old 0 -> new 1  (spatial sigma^0 -> spatial sigma^1)
    old 1 -> new 2  (spatial sigma^1 -> spatial sigma^2)
    old 2 -> new 3  (spatial sigma^2 -> spatial sigma^3)
    old 3 -> new 0  (fiber/L dtau   -> time N(t) dt)

--- ADAPTER BOUNDARY POLICY ---

This module is the ONLY permitted entry point for converting old-DPPU
(Euclidean) topology data into DPPU-LZ (Lorentzian) native convention.

Permitted import locations:
  - dppu/topology/ modules
  - topology-layer check scripts (Phase 1A, 1C, ...)

FORBIDDEN import locations:
  - dppu/engine/   (metric, levi_civita, ec_connection, contortion, pipeline)
  - dppu/curvature/ (riemann, weyl, hodge, pontryagin, ...)
  - dppu/action/
  - dppu/torsion/

The core pipeline speaks DPPU-LZ (indices 0,1,2,3 with 0=time) natively.
Legacy conversion belongs at the topology/data-ingestion boundary only.

--- SIGN NOTE ---

The index permutation sign sgn(f) = -1 (Phase 1A result) applies to
orientation-sensitive objects (epsilon, Hodge dual).  Structure constant
VALUES are preserved by this adapter; sign is NOT re-introduced here.
"""

from typing import Dict, Tuple, Any


# ---------------------------------------------------------------------------
# Index mapping helpers (thin wrappers around legacy_index module)
# ---------------------------------------------------------------------------

def old_to_lz_index(a_old: int) -> int:
    """Map single frame index from old-DPPU to LZ-native: a_new = (a_old+1)%4."""
    from ..engine.legacy_index import old_to_lz_index as _f
    return _f(a_old)


def lz_to_old_index(a_lz: int) -> int:
    """Map single frame index from LZ-native to old-DPPU: a_old = (a_lz+3)%4."""
    from ..engine.legacy_index import lz_to_old_index as _f
    return _f(a_lz)


# ---------------------------------------------------------------------------
# Structure constant remapping
# ---------------------------------------------------------------------------

def remap_structure_constants_to_lz(
    C_old: Dict[Tuple[int, int, int], Any],
    dim: int = 4,
) -> Dict[Tuple[int, int, int], Any]:
    """
    Remap a sparse structure-constant dict from old-DPPU to LZ-native convention.

    Args:
        C_old : dict  {(a, b, c): value}  in old-DPPU (Euclidean) index convention.
                Only nonzero entries need to be present.
        dim   : frame dimension (default 4).

    Returns:
        C_lz  : dict  {(a_new, b_new, c_new): value}  in LZ-native index convention.
                Values are unchanged; only keys are permuted.

    Notes:
        - C^a_{bc} antisymmetry in (b,c) is preserved by the permutation.
        - Values are copied as-is; no sign factor is applied here.
          (sgn(f)=-1 is an orientation sign for epsilon/Hodge, not for C itself.)
    """
    C_lz: Dict[Tuple[int, int, int], Any] = {}
    for (a_old, b_old, c_old), val in C_old.items():
        a_new = old_to_lz_index(a_old)
        b_new = old_to_lz_index(b_old)
        c_new = old_to_lz_index(c_old)
        C_lz[(a_new, b_new, c_new)] = val
    return C_lz


# ---------------------------------------------------------------------------
# Per-topology LZ-native structure constant data
# ---------------------------------------------------------------------------

def get_t3_lz_structure_constants(dim: int = 4) -> Dict[Tuple[int, int, int], Any]:
    """
    Return LZ-native structure constants for T3 (all zero).

    T3 is a flat manifold: C^a_{bc} = 0 for all indices in both conventions.
    This function returns an empty dict (no nonzero entries).
    """
    return {}


def get_s3_lz_index_mapping() -> Dict[str, Any]:
    """
    Return the S3 old->LZ index-block mapping for documentation and auditing.

    Old-DPPU S3:
        spatial block  a,b,c in {0,1,2}
        fiber          a = 3

    LZ-native S3:
        time           a = 0
        spatial block  a,b,c in {1,2,3}

    The structure constant VALUES are preserved; only the key indices change.
    """
    return {
        "old_spatial_block": {0, 1, 2},
        "new_spatial_block": {1, 2, 3},
        "old_fiber_index":   3,
        "new_time_index":    0,
        "description": (
            "old spatial {0,1,2} -> new spatial {1,2,3}; "
            "old fiber 3 -> new time 0"
        ),
    }


def get_nil3_lz_key_mapping() -> Dict[Tuple[int, int, int], Tuple[int, int, int]]:
    """
    Return the Nil3 nonzero-key mapping from old-DPPU to LZ-native convention.

    Old-DPPU nonzero keys (spatial block {0,1,2}, antisymmetric):
        (2, 0, 1)  <->  (2, 1, 0)

    LZ-native nonzero keys (spatial block {1,2,3}):
        (3, 1, 2)  <->  (3, 2, 1)

    Values are unchanged.
    """
    return {
        (2, 0, 1): (3, 1, 2),
        (2, 1, 0): (3, 2, 1),
    }


def verify_t3_lz_all_zero(
    C_lz: Dict[Tuple[int, int, int], Any]
) -> bool:
    """
    Verify that all structure constants for T3 are zero in LZ-native convention.

    Args:
        C_lz: dict {(a,b,c): value} or MutableDenseNDimArray; nonzero entries only.

    Returns:
        True if all values are zero (or dict is empty).
    """
    from sympy import S
    if hasattr(C_lz, '__iter__') and not hasattr(C_lz, '__getitem__'):
        # MutableDenseNDimArray or similar
        return all(v == S.Zero for v in C_lz)
    if isinstance(C_lz, dict):
        return all(v == S.Zero or v == 0 for v in C_lz.values())
    # MutableDenseNDimArray: iterate over all elements
    return all(v == S.Zero or v == 0 for v in C_lz)


def verify_s3_spatial_keys_in_lz_block(
    C_lz: Dict[Tuple[int, int, int], Any]
) -> Tuple[bool, list]:
    """
    Verify that all nonzero S3 structure constants have keys only in the LZ spatial block {1,2,3}.

    Returns:
        (passed, violations)
        violations: list of (a, b, c) keys that fall outside the {1,2,3} block.
    """
    violations = []
    for (a, b, c), val in C_lz.items():
        if val == 0:
            continue
        # All spatial indices should be in {1,2,3}; no index should be 0 (time)
        if 0 in (a, b, c):
            violations.append((a, b, c))
    return len(violations) == 0, violations


def verify_nil3_lz_keys(
    C_lz: Dict[Tuple[int, int, int], Any]
) -> Tuple[bool, list]:
    """
    Verify that Nil3 nonzero structure constants use the expected LZ-native keys.

    Expected nonzero keys: (3,1,2) and (3,2,1).

    Returns:
        (passed, unexpected_keys)
    """
    expected = {(3, 1, 2), (3, 2, 1)}
    unexpected = []
    for (a, b, c), val in C_lz.items():
        if val == 0:
            continue
        if (a, b, c) not in expected:
            unexpected.append((a, b, c))
    return len(unexpected) == 0, unexpected


# ---------------------------------------------------------------------------
# Dense array remap (for use in base_topology.py after _build_structure_constants)
# ---------------------------------------------------------------------------

def remap_dense_structure_constants_to_lz(C_old, dim: int = 4):
    """
    Remap a dense MutableDenseNDimArray C^a_{bc} from old-DPPU convention
    to LZ-native convention.

    Index permutation: a_new = (a_old + 1) % 4
    Values are preserved; no orientation sign is applied.

    This helper is called by base_topology.step_E4_2_metric_and_frame
    after _build_structure_constants when signature="lorentzian" and
    frame_convention="lz_native".  It MUST NOT be called from engine/
    curvature/action/torsion core (adapter boundary policy).

    Args:
        C_old: MutableDenseNDimArray of shape (dim, dim, dim)
        dim:   Frame dimension (default 4)

    Returns:
        C_new: MutableDenseNDimArray with keys permuted to LZ-native
    """
    from sympy import S as _S
    from sympy.tensor.array import MutableDenseNDimArray

    C_new = MutableDenseNDimArray.zeros(dim, dim, dim)
    for a_old in range(dim):
        for b_old in range(dim):
            for c_old in range(dim):
                val = C_old[a_old, b_old, c_old]
                if val != _S.Zero and val != 0:
                    a_new = old_to_lz_index(a_old)
                    b_new = old_to_lz_index(b_old)
                    c_new = old_to_lz_index(c_old)
                    C_new[a_new, b_new, c_new] = val
    return C_new
