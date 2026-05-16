"""
Legacy Index Adapter: Euclidean DPPU -> Lorentzian DPPU-LZ
==========================================================

Maps frame indices between the two conventions:

  Old DPPU (Euclidean):
    a=0,1,2  spatial frame (SU(2) structure constants)
    a=3      fiber direction (Euclidean S1, length L)
    metric:  diag(+1,+1,+1,+1)

  DPPU-LZ (Lorentzian):
    a=0      time direction (Lorentzian, lapse N(t))
    a=1,2,3  spatial frame (SU(2) structure constants)
    metric:  diag(-1,+1,+1,+1)

Permutation map:
    old_to_lz: a_new = (a_old + 1) % 4
    lz_to_old: a_old = (a_new + 3) % 4

Map table:
    old 0 -> new 1  (1st spatial: sigma^0 -> sigma^1)
    old 1 -> new 2  (2nd spatial: sigma^1 -> sigma^2)
    old 2 -> new 3  (3rd spatial: sigma^2 -> sigma^3)
    old 3 -> new 0  (fiber/time:  L dtau  -> N(t) dt)

--- ADAPTER BOUNDARY POLICY ---

This module is the ONLY permitted entry point for converting old-DPPU
(Euclidean) data into DPPU-LZ (Lorentzian) convention.

Permitted usage locations:
  - dppu/topology/ layer, when loading old-DPPU topology data
  - scripts that read legacy-format checkpoint/validation files
  - Phase 1A check scripts

FORBIDDEN usage:
  - Inside dppu/engine/ core pipeline modules (ec_connection, levi_civita,
    contortion, metric, pipeline)
  - Inside dppu/curvature/ modules (hodge, riemann, pontryagin, weyl, etc.)
  - Inside dppu/action/ modules (ec_action, lagrangian, etc.)
  - Inside dppu/torsion/ modules

The core pipeline must speak DPPU-LZ (0,1,2,3) natively.
Legacy conversion belongs at the topology/data-ingestion boundary only.

--- ORIENTATION-SENSITIVE OBJECTS ---

IMPORTANT: This adapter handles two distinct classes of objects:

  1. Ordinary tensor components:
     Index relabeling only.  No extra sign.
     Use: permute_component_key(), permute_tensor_components()

  2. Orientation-sensitive objects (epsilon symbol, Hodge dual, etc.):
     The permutation f: (0,1,2,3) -> (1,2,3,0) is a 4-cycle (odd permutation).
     Its sign is sgn(f) = (-1)^(4-1) = -1.
     When comparing old-convention epsilon with new-convention epsilon,
     a factor of sgn(f) = -1 must be applied explicitly.
     Use: legacy_epsilon_permutation_sign()

  DO NOT use permute_tensor_components() blindly for epsilon symbols.

--- L vs N(t) ---

The old fiber coordinate e^3_old = L dtau maps to e^0_new = N(t) dt.
L (Euclidean S1 circumference, constant) and N(t) (Lorentzian lapse,
time-dependent) are physically distinct objects.  They must NOT be
identified in code.  This adapter records the index correspondence only.

Reference: DPPU-LZ_legacy_index_map_v1.1.md, DPPU-LZ_conventions_v1.1.md
"""

from typing import Dict, Tuple


# ---------------------------------------------------------------------------
# Core index maps
# ---------------------------------------------------------------------------

def validate_frame_index(a: int) -> None:
    """Raise ValueError unless a is one of 0, 1, 2, 3."""
    if a not in (0, 1, 2, 3):
        raise ValueError(
            f"Frame index {a!r} is out of range; expected one of {{0, 1, 2, 3}}."
        )


def old_to_lz_index(a_old: int) -> int:
    """Map Euclidean DPPU legacy frame index to Lorentzian DPPU-LZ frame index.

    Permutation: a_new = (a_old + 1) % 4

    Raises:
        ValueError: If a_old is not in {0, 1, 2, 3}.

    Examples:
        >>> old_to_lz_index(3)  # fiber -> time
        0
        >>> old_to_lz_index(0)  # 1st spatial
        1
    """
    validate_frame_index(a_old)
    return (a_old + 1) % 4


def lz_to_old_index(a_lz: int) -> int:
    """Map Lorentzian DPPU-LZ frame index back to Euclidean DPPU legacy frame index.

    Inverse permutation: a_old = (a_lz + 3) % 4

    Raises:
        ValueError: If a_lz is not in {0, 1, 2, 3}.

    Examples:
        >>> lz_to_old_index(0)  # time -> fiber
        3
        >>> lz_to_old_index(1)  # 1st spatial
        0
    """
    validate_frame_index(a_lz)
    return (a_lz + 3) % 4


# ---------------------------------------------------------------------------
# Tuple helpers
# ---------------------------------------------------------------------------

def old_to_lz_indices(indices: Tuple[int, ...]) -> Tuple[int, ...]:
    """Map a tuple of legacy indices to LZ-native indices."""
    return tuple(old_to_lz_index(a) for a in indices)


def lz_to_old_indices(indices: Tuple[int, ...]) -> Tuple[int, ...]:
    """Map a tuple of LZ-native indices to legacy indices."""
    return tuple(lz_to_old_index(a) for a in indices)


# ---------------------------------------------------------------------------
# Permutation sign
# ---------------------------------------------------------------------------

def legacy_permutation_sign_4d() -> int:
    """Return the orientation sign of the old->LZ 4D permutation.

    The map f: (0,1,2,3) -> (1,2,3,0) is a single 4-cycle.
    A k-cycle has sign (-1)^(k-1); for k=4: sgn(f) = (-1)^3 = -1.

    Verification by inversion count:
      Sequence [1,2,3,0] has inversions (1>0),(2>0),(3>0) -> 3 inversions.
      Sign = (-1)^3 = -1.

    Returns: -1
    """
    perm = [old_to_lz_index(a) for a in range(4)]  # [1, 2, 3, 0]
    inversions = sum(
        1
        for i in range(4)
        for j in range(i + 1, 4)
        if perm[i] > perm[j]
    )
    return 1 if inversions % 2 == 0 else -1


# ---------------------------------------------------------------------------
# Epsilon orientation sign (orientation-sensitive object)
# ---------------------------------------------------------------------------

def legacy_epsilon_permutation_sign() -> int:
    """Return sign needed when comparing canonical epsilon symbols under
    the old->LZ permutation.

    old canonical order: (0,1,2,3)
    mapped order:        (1,2,3,0)
    sign:                -1

    Usage:
        epsilon_symbol_new(1,2,3,0) = sgn(f) * epsilon_symbol_old(0,1,2,3)
                                    = -1 * 1 = -1

    DO NOT confuse with ordinary tensor index relabeling.
    This sign must be applied explicitly when converting epsilon-dependent
    expressions (Hodge dual, Pontryagin density, NY term) from old to new
    convention.

    Returns: -1
    """
    return legacy_permutation_sign_4d()


# ---------------------------------------------------------------------------
# Component dictionary helpers (ordinary tensors only)
# ---------------------------------------------------------------------------

def permute_component_key(
    key: Tuple[int, ...],
    direction: str = "old_to_lz",
) -> Tuple[int, ...]:
    """Permute a tensor component key between legacy and LZ conventions.

    This function only relabels indices and does NOT apply orientation signs.
    For epsilon symbols or other orientation-sensitive objects, use
    legacy_epsilon_permutation_sign() to obtain the additional sign factor.

    Args:
        key: Tuple of frame indices (each in {0, 1, 2, 3}).
        direction: "old_to_lz" (default) or "lz_to_old".

    Returns:
        Remapped key tuple.

    Examples:
        >>> permute_component_key((0, 1, 2, 3), "old_to_lz")
        (1, 2, 3, 0)
        >>> permute_component_key((1, 2, 3, 0), "lz_to_old")
        (0, 1, 2, 3)
    """
    if direction == "old_to_lz":
        fn = old_to_lz_index
    elif direction == "lz_to_old":
        fn = lz_to_old_index
    else:
        raise ValueError(
            f"Unknown direction {direction!r}; expected 'old_to_lz' or 'lz_to_old'."
        )
    return tuple(fn(a) for a in key)


def permute_tensor_components(
    components: Dict[Tuple[int, ...], object],
    direction: str = "old_to_lz",
) -> Dict[Tuple[int, ...], object]:
    """Relabel component dictionary keys between legacy and LZ conventions.

    This is for ORDINARY tensor components only.
    Do NOT use this blindly for epsilon symbols or other orientation-sensitive
    objects; those require an additional sgn(f) = -1 factor via
    legacy_epsilon_permutation_sign().

    Args:
        components: Dict mapping index tuples to component values.
        direction: "old_to_lz" (default) or "lz_to_old".

    Returns:
        New dict with relabeled keys.  Values are unchanged.
    """
    return {
        permute_component_key(key, direction): val
        for key, val in components.items()
    }
