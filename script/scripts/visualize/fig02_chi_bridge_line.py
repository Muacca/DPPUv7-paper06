"""Generate Fig 2: Diagnostic bridge — C_topology = -9*chi line.

Data from chi_bridge_symbolic_20260513_230305.log (all PASS).
The line C = -9*chi is the family-level identity; the 4 points are specializations.

Output: LaTeX/figures/fig02_chi_bridge_line.png
"""

from __future__ import annotations

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[2]
FIGURE_DIR = PROJECT_ROOT / "LaTeX" / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# ── Verified topology data (from chi_bridge_symbolic log, all PASS) ──────────
TOPOLOGY_DATA = [
    #  label      chi      C       marker  color       label_ax_x  label_ax_y
    ("$T^3$",    0,       0,      'o', '#222222',    0.42,        0.78),
    ("$Nil^3$",  -1/12,   3/4,   '^', '#1a6fc4',    0.42,        0.88),
    ("$Sol^3$",  -1/3,    3,     's', '#c46a1a',    0.42,        0.97),
    ("$S^3$",    4,      -36,   'D', '#b01818',    0.88,        0.16),
]

XLO, XHI = -0.55, 4.8
YLO, YHI = -42, 8


def _frac(x: float) -> str:
    f = Fraction(x).limit_denominator(24)
    if f.denominator == 1:
        return str(f.numerator)
    return f'\\frac{{{f.numerator}}}{{{f.denominator}}}'


def main() -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5.2))

    # ── Reference line C = -9*chi ────────────────────────────────────────────
    chi_line = np.array([XLO, XHI])
    ax.plot(chi_line, -9 * chi_line,
            color='#555555', lw=1.2, linestyle='-', zorder=1)

    # ── Data points ──────────────────────────────────────────────────────────
    for label, chi, C, marker, color, ax_x, ax_y in TOPOLOGY_DATA:
        ax.scatter(chi, C, marker=marker, color=color, s=75, zorder=5)
        ax.annotate(
            f'{label}\n$({_frac(chi)},\\ {_frac(C)})$',
            xy=(chi, C),                  # data coordinates
            xytext=(ax_x, ax_y),          # axes fraction
            xycoords='data',
            textcoords='axes fraction',
            fontsize=8.5, color=color, va='center',
            arrowprops=dict(
                arrowstyle='-',
                color=color, lw=0.8,
                connectionstyle='arc3,rad=0.15',
            ),
        )

    # ── Axes ─────────────────────────────────────────────────────────────────
    ax.axhline(0, color='#aaaaaa', lw=0.7, zorder=0)
    ax.axvline(0, color='#aaaaaa', lw=0.7, zorder=0)
    ax.grid(True, which='major', color='#dddddd', lw=0.6, zorder=0)
    ax.grid(True, which='minor', color='#eeeeee', lw=0.4, zorder=0)

    ax.set_xlim(XLO, XHI)
    ax.set_ylim(YLO, YHI)
    ax.set_xlabel(r'$\chi$', fontsize=11)
    ax.set_ylabel(r'$C_{\rm topology}$', fontsize=11)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.25))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(9))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(3))

    # ── Legend ───────────────────────────────────────────────────────────────
    line_handle = plt.Line2D([0], [0], color='#555555', lw=1.2,
                             label=r'$C + 9\chi = 0$  (family-level identity)')
    marker_handles = [
        plt.Line2D([0], [0], marker=m, color='none', markerfacecolor=c,
                   markersize=7, label=lbl)
        for lbl, _chi, _C, m, c, *_ in TOPOLOGY_DATA
    ]
    ax.legend(handles=[line_handle] + marker_handles,
              fontsize=8.5, loc='lower left', framealpha=0.92,
              edgecolor='#cccccc')

    fig.tight_layout()
    out = FIGURE_DIR / "fig02_chi_bridge_line.png"
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
