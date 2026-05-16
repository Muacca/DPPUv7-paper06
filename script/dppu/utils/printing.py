"""
DPPU Printing Utilities
=======================

Console output helpers for consistent formatting across scripts.

These functions are intentionally simple so that standalone scripts
can reproduce the same look without importing this module.

Author: Muacca
"""


def hline(char: str = '=', width: int = 70) -> None:
    """Print a horizontal rule."""
    print(char * width)


def print_header(title: str, width: int = 70) -> None:
    """Print a section header with surrounding rules."""
    print()
    hline(width=width)
    print(f"  {title}")
    hline(width=width)


def print_sub(title: str) -> None:
    """Print a sub-section title."""
    print()
    print(f"  \u2500\u2500 {title} \u2500\u2500")
    print()
