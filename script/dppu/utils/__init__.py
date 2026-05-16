"""
DPPU Utils: Shared Utility Functions
====================================

Provides common utility functions used across the package.

Modules:
- epsilon: Totally antisymmetric epsilon symbols (3D and 4D)
- symbolic: Symbolic computation helpers (prove_zero, witness search)
- logger: ComputationLogger and NullLogger
- printing: Console output helpers (hline, print_header, print_sub)
- tee_logger: Stdout tee for script log capture
- io: CSV read/write helpers
"""

from .epsilon import epsilon_3d, epsilon_4d, epsilon_nd
from .symbolic import (
    prove_zero,
    find_nonzero_witness,
    generate_test_points,
    derivative_inventory,
    nonzero_component_count,
)
from .logger import ComputationLogger, NullLogger
from .printing import hline, print_header, print_sub
from .io import (
    now_utc_iso,
    read_text,
    write_text,
    read_json,
    write_json,
    read_csv,
    write_csv,
    upsert_csv_rows,
)

__all__ = [
    'epsilon_3d',
    'epsilon_4d',
    'epsilon_nd',
    'prove_zero',
    'find_nonzero_witness',
    'generate_test_points',
    'derivative_inventory',
    'nonzero_component_count',
    'ComputationLogger',
    'NullLogger',
    'hline',
    'print_header',
    'print_sub',
    'now_utc_iso',
    'read_text',
    'write_text',
    'read_json',
    'write_json',
    'read_csv',
    'write_csv',
    'upsert_csv_rows',
]
