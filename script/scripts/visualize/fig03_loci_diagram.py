"""Generate Fig 3: Distinguished loci in the five-parameter family.

Two-panel schematic:
  Left  — class-A (a, c) plane  [b = a for all four topologies]
  Right — class-B (u, v) plane  [c fixed per topology]

Data: verified loci from chi_bridge_symbolic log (all PASS).

Output: LaTeX/figures/fig03_loci_diagram.png
"""

from __future__ import annotations

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[2]
FIGURE_DIR = PROJECT_ROOT / "LaTeX" / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

# ── Topology parameter data (a, b, c, u, v) ──────────────────────────────────
# All from Table 5b; b = a for every topology.
TOPOLOGIES = {
    "$T^3$":   dict(a=0,  b=0,  c=0,  u=0,  v=0,  chi=0,      marker='o', color='#222222'),
    "$Nil^3$": dict(a=0,  b=0,  c=-1, u=0,  v=0,  chi=-1/12,  marker='^', color='#1a6fc4'),
    "$S^3$":   dict(a=4,  b=4,  c=4,  u=0,  v=0,  chi=4,      marker='D', color='#b01818'),
    "$Sol^3$": dict(a=0,  b=0,  c=0,  u=1,  v=-1, chi=-1/3,   marker='s', color='#c46a1a'),
}


def _label_chi(chi: float) -> str:
    from fractions import Fraction
    f = Fraction(chi).limit_denominator(24)
    if f.denominator == 1:
        s = str(f.numerator)
    else:
        s = f'\\frac{{{f.numerator}}}{{{f.denominator}}}'
    return f'$\\chi={s}$'


def panel_A(ax: plt.Axes) -> None:
    """Left panel: class-A (a, c) space. b = a for all topologies."""
    # Isotropic diagonal a = b = c → in (a, c) plane: c = a
    a_diag = np.array([-0.3, 4.8])
    ax.plot(a_diag, a_diag, linestyle='--', color='#888888', lw=1.2,
            label='$a = b = c$ (isotropic diagonal)', zorder=1)

    # Axes guides
    ax.axhline(0, color='#cccccc', lw=0.7, zorder=0)
    ax.axvline(0, color='#cccccc', lw=0.7, zorder=0)

    offset = {
        "$T^3$":   (0.1,  1.0),
        "$Nil^3$": (0.5,  -0.1),   # to the left of the axis
        "$S^3$":   (-0.65,  0.40),
        "$Sol^3$": (0.4, 0.0),    # below T^3 in class-A
    }
    plotted_origin = False
    for name, p in TOPOLOGIES.items():
        x, y = p['a'], p['c']
        # T^3 and Sol^3 both at origin in class-A — plot once, label both
        if x == 0 and y == 0:
            if not plotted_origin:
                ax.scatter(x, y, marker=p['marker'], color=p['color'],
                           s=70, zorder=5)
                plotted_origin = True
            else:
                # Sol^3 at origin in class-A: use its own offset entry
                dx, dy = offset[name]
                ax.annotate(
                    f'{name}\n{_label_chi(p["chi"])}',
                    xy=(x, y), xytext=(x + dx, y + dy),
                    fontsize=8, color=p['color'], va='top',
                    arrowprops=dict(arrowstyle='-', color=p['color'], lw=0.7),
                )
                continue
        else:
            ax.scatter(x, y, marker=p['marker'], color=p['color'],
                       s=70, zorder=5)

        dx, dy = offset[name]
        ax.annotate(
            f'{name}\n{_label_chi(p["chi"])}',
            xy=(x, y), xytext=(x + dx, y + dy),
            fontsize=8, color=p['color'], va='center',
            arrowprops=dict(arrowstyle='-', color=p['color'], lw=0.7),
        )

    ax.set_xlim(-0.5, 5.2)
    ax.set_ylim(-1.7, 5.0)
    ax.set_xlabel('$a$ (= $b$)', fontsize=10)
    ax.set_ylabel('$c$', fontsize=10)
    ax.set_title('(a)  class-A subspace $(a, b, c)$', fontsize=10, pad=6)

    handles = [
        plt.Line2D([0], [0], linestyle='--', color='#888888', lw=1.2,
                   label='$a = b = c$ isotropic diagonal'),
    ] + [
        plt.Line2D([0], [0], marker=p['marker'], color='none',
                   markerfacecolor=p['color'], markersize=7, label=name)
        for name, p in TOPOLOGIES.items()
    ]
    ax.legend(handles=handles, fontsize=7.5, loc='upper left',
              framealpha=0.9, edgecolor='#cccccc')
    ax.grid(True, color='#eeeeee', lw=0.5, zorder=0)


def panel_B(ax: plt.Axes) -> None:
    """Right panel: class-B (u, v) space."""
    # u + v = 0 line → v = -u
    u_line = np.array([-0.3, 1.5])
    ax.plot(u_line, -u_line, linestyle=':', color='#888888', lw=1.2,
            label='$u + v = 0$', zorder=1)

    ax.axhline(0, color='#cccccc', lw=0.7, zorder=0)
    ax.axvline(0, color='#cccccc', lw=0.7, zorder=0)

    # T^3, Nil^3, S^3 all have u=v=0 → cluster at origin
    origin_labels = [
        (name, p) for name, p in TOPOLOGIES.items() if p['u'] == 0 and p['v'] == 0
    ]
    # Plot a single point for the origin cluster
    ax.scatter(0, 0, marker='o', color='#222222', s=70, zorder=5)

    # Staggered labels for the three origin topologies
    offsets_origin = [(0.0, -0.14), (0.0, -0.30), (0.0, -0.46)]
    for (name, p), (dx, dy) in zip(origin_labels, offsets_origin):
        ax.text(dx, dy, f'{name}  $\\chi={_chi_str(p["chi"])}$',
                ha='right', va='center', fontsize=8, color=p['color'])

    # Sol^3 at (u=1, v=-1)
    p_sol = TOPOLOGIES["$Sol^3$"]
    ax.scatter(p_sol['u'], p_sol['v'], marker=p_sol['marker'],
               color=p_sol['color'], s=70, zorder=5)
    ax.annotate(
        f'$Sol^3$\n$\\chi={_chi_str(p_sol["chi"])}$\n$(u,v)=(1,-1)$',
        xy=(p_sol['u'], p_sol['v']),
        xytext=(p_sol['u'] + 0.12, p_sol['v'] - 0.18),
        fontsize=8, color=p_sol['color'], va='top',
        arrowprops=dict(arrowstyle='-', color=p_sol['color'], lw=0.7),
    )

    ax.set_xlim(-0.55, 1.65)
    ax.set_ylim(-1.55, 0.55)
    ax.set_xlabel('$u$', fontsize=10)
    ax.set_ylabel('$v$', fontsize=10)
    ax.set_title('(b)  class-B subspace $(u, v)$', fontsize=10, pad=6)

    ax.legend(fontsize=7.5, loc='upper right', framealpha=0.9, edgecolor='#cccccc')
    ax.grid(True, color='#eeeeee', lw=0.5, zorder=0)


def _chi_str(chi: float) -> str:
    from fractions import Fraction
    f = Fraction(chi).limit_denominator(24)
    if f.denominator == 1:
        return str(f.numerator)
    return f'\\frac{{{f.numerator}}}{{{f.denominator}}}'


def main() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    fig.subplots_adjust(wspace=0.35)

    panel_A(axes[0])
    panel_B(axes[1])

    out = FIGURE_DIR / "fig03_loci_diagram.png"
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
