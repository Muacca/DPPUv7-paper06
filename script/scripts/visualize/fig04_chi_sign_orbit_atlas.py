"""Generate Fig 4: chi-sign reduced orbit atlas schematic.

Three-panel qualitative q–t diagram:
  (a) chi > 0  (S^3, chi=4):   no real trajectory — hatched region
  (b) chi = 0  (T^3, chi=0):   static / degenerate — horizontal line
  (c) chi < 0  (Nil^3/Sol^3): monotonic expansion + singular-approach collapse

Slopes in panel (c):
  Nil^3: qdot = ±sqrt(-chi) = ±sqrt(1/12) = ±1/(2*sqrt(3)) ≈ ±0.2887
  Sol^3: qdot = ±sqrt(-chi) = ±sqrt(1/3)  = ±1/sqrt(3)     ≈ ±0.5774

Output: LaTeX/figures/fig04_chi_sign_orbit_atlas.png
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

T_MAX = 3.0
Q0 = 1.0          # initial q at t = 0
T = np.linspace(0, T_MAX, 400)


def panel_a(ax: plt.Axes) -> None:
    """chi > 0: S^3. No real trajectory — qdot^2 = -chi < 0."""
    ax.set_xlim(0, T_MAX)
    ax.set_ylim(0, 3.0)

    # Hatched region filling the entire panel
    ax.fill_between(T, 0, 3.0, color='#cccccc', alpha=0.55,
                    hatch='////', linewidth=0.0)
    ax.text(T_MAX / 2, 1.5,
            'no real trajectory\n$\\dot{q}^2 = -\\chi < 0$',
            ha='center', va='center', fontsize=9,
            color='#444444',
            bbox=dict(boxstyle='round,pad=0.35', facecolor='white',
                      edgecolor='#888888', alpha=0.85))

    ax.set_xlabel('$t$', fontsize=10)
    ax.set_ylabel('$q$', fontsize=10)
    ax.set_title(r'(a)  $\chi > 0$ — $S^3\,(\chi=4)$', fontsize=10, pad=6)
    ax.tick_params(labelleft=False, left=False)


def panel_b(ax: plt.Axes) -> None:
    """chi = 0: T^3. qdot = 0 → q = const (static / degenerate)."""
    ax.set_xlim(0, T_MAX)
    ax.set_ylim(0, 3.0)

    q_static = Q0 + 1.0    # arbitrary constant height for visual clarity
    ax.plot(T, np.full_like(T, q_static), color='#222222', lw=2.0,
            label='$q(t) = \\mathrm{const}$')

    ax.annotate('static / degenerate\n$\\dot{q} = 0$',
                xy=(T_MAX * 0.6, q_static),
                xytext=(T_MAX * 0.6, q_static + 0.55),
                fontsize=8.5, ha='center', va='bottom',
                arrowprops=dict(arrowstyle='->', color='#222222', lw=0.9))

    ax.set_xlabel('$t$', fontsize=10)
    ax.set_ylabel('$q$', fontsize=10)
    ax.set_title(r'(b)  $\chi = 0$ — $T^3\,(\chi=0)$', fontsize=10, pad=6)
    ax.tick_params(labelleft=False, left=False)


def panel_c(ax: plt.Axes) -> None:
    """chi < 0: Nil^3 and Sol^3. Two sheets each."""
    ax.set_xlim(0, T_MAX)
    ax.set_ylim(0, 3.5)

    chi_nil = -1 / 12
    chi_sol = -1 / 3
    slope_nil = np.sqrt(-chi_nil)    # ≈ 0.2887
    slope_sol = np.sqrt(-chi_sol)    # ≈ 0.5774

    # Reference point: q(0) = 1.5 (mid-panel)
    q0 = 1.5

    # Nil^3 expanding (slope > 0): dashed
    ax.plot(T, q0 + slope_nil * T, color='#1a6fc4', lw=1.8,
            linestyle='--', label=f'$Nil^3$ expanding  $\\dot{{q}}_0={slope_nil:.3f}$')
    # Nil^3 collapsing (slope < 0): dashed — clip at Q_FLOOR
    q_coll_nil = q0 - slope_nil * T
    mask_nil = q_coll_nil > 0.05
    ax.plot(T[mask_nil], q_coll_nil[mask_nil], color='#1a6fc4', lw=1.8,
            linestyle='--', label=f'$Nil^3$ collapse  $\\dot{{q}}_0={-slope_nil:.3f}$')

    # Sol^3 expanding: solid
    ax.plot(T, q0 + slope_sol * T, color='#c46a1a', lw=1.8,
            linestyle='-', label=f'$Sol^3$ expanding  $\\dot{{q}}_0={slope_sol:.3f}$')
    # Sol^3 collapsing: solid
    q_coll_sol = q0 - slope_sol * T
    mask_sol = q_coll_sol > 0.05
    ax.plot(T[mask_sol], q_coll_sol[mask_sol], color='#c46a1a', lw=1.8,
            linestyle='-', label=f'$Sol^3$ collapse  $\\dot{{q}}_0={-slope_sol:.3f}$')

    # Arrows indicating time direction (expanding sheets)
    for slope, color, ls in [
        (slope_nil, '#1a6fc4', '--'),
        (slope_sol, '#c46a1a', '-'),
    ]:
        t_arr = T_MAX * 0.55
        q_arr = q0 + slope * t_arr
        ax.annotate('',
                    xy=(t_arr + 0.25, q0 + slope * (t_arr + 0.25)),
                    xytext=(t_arr, q_arr),
                    arrowprops=dict(arrowstyle='->', color=color, lw=1.2))

    # Label "expanding" / "collapsing" inside the panel
    ax.text(T_MAX * 0.50, q0 + slope_sol * T_MAX * 0.60 + 0.12,
            'expanding', ha='center', va='bottom', fontsize=8, color='#555555')
    ax.text(T_MAX * 0.40, q0 - slope_sol * T_MAX * 0.42 - 0.12,
            'collapse\n($q \\to 0$)', ha='center', va='top', fontsize=8, color='#555555')

    ax.set_xlabel('$t$', fontsize=10)
    ax.set_ylabel('$q$', fontsize=10)
    ax.set_title(r'(c)  $\chi < 0$ — $Nil^3,\,Sol^3$', fontsize=10, pad=6)
    ax.tick_params(labelleft=False, left=False)

    legend_handles = [
        plt.Line2D([0], [0], color='#1a6fc4', lw=1.8, linestyle='--',
                   label=r'$Nil^3$  ($\chi = -1/12$,  $|\dot{q}| \approx 0.289$)'),
        plt.Line2D([0], [0], color='#c46a1a', lw=1.8, linestyle='-',
                   label=r'$Sol^3$  ($\chi = -1/3$,   $|\dot{q}| \approx 0.577$)'),
    ]
    ax.legend(handles=legend_handles, fontsize=8.0, loc='upper left',
              framealpha=0.92, edgecolor='#cccccc')


def main() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
    fig.subplots_adjust(wspace=0.30)

    panel_a(axes[0])
    panel_b(axes[1])
    panel_c(axes[2])

    # Shared y-label note
    for ax in axes:
        ax.tick_params(axis='x', labelsize=8.5)

    out = FIGURE_DIR / "fig04_chi_sign_orbit_atlas.png"
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
