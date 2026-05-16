"""Generate Fig 1: Three-level architecture of chi-controlled topology universality.

Output: LaTeX/figures/fig01_three_level_architecture.png
"""

from __future__ import annotations

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT_ROOT = HERE.parents[2]   # .../20_DPPUv7-paper06/
FIGURE_DIR = PROJECT_ROOT / "LaTeX" / "figures"
FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def draw_box(ax, x, y, text, fc, ec, fontsize=9.5, lw=1.4):
    return ax.text(
        x, y, text,
        ha='center', va='center', fontsize=fontsize,
        bbox=dict(boxstyle='round,pad=0.45', facecolor=fc, edgecolor=ec, lw=lw),
        multialignment='center',
        linespacing=1.5,
    )


def draw_arrow(ax, x1, y1, x2, y2, color='#2c5fa0', lw=1.4, linestyle='solid'):
    ax.annotate(
        '', xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle='->', color=color, lw=lw,
            linestyle=linestyle,
            mutation_scale=12,
        ),
    )


def main() -> None:
    fig, ax = plt.subplots(figsize=(8, 6.2))
    ax.axis('off')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)

    # ── Top box: Five-parameter family ──────────────────────────────────────
    draw_box(
        ax, 5, 9.05,
        "Five-parameter family $(a,\\,b,\\,c,\\,u,\\,v)$\n"
        "$S^3\\,(\\chi{=}4),\\ T^3\\,(\\chi{=}0),\\ "
        "Nil^3\\,(\\chi{=}{-}\\frac{1}{12}),\\ Sol^3\\,(\\chi{=}{-}\\frac{1}{3})$",
        fc='#d6e8f5', ec='#2c5fa0', fontsize=9.8,
    )

    # Arrow: top box → chi box
    draw_arrow(ax, 5, 8.42, 5, 7.78, color='#2c5fa0')
    ax.text(5.18, 8.10, 'curvature reduction\n'
            r'$\chi = {}^{(3)}R\,q^2/6$',
            ha='left', va='center', fontsize=8, color='#2c5fa0', fontstyle='italic',
            linespacing=1.4)

    # ── Middle box: Scalar χ ────────────────────────────────────────────────
    draw_box(
        ax, 5, 7.30,
        r"Scalar $\chi$",
        fc='#fff6c8', ec='#b07800', fontsize=11, lw=1.6,
    )

    # Vertical stem from χ box to junction bar
    ax.plot([5, 5], [6.85, 6.40], color='#2c5fa0', lw=1.4, solid_capstyle='butt')

    # Horizontal junction bar
    ax.plot([1.9, 8.1], [6.40, 6.40], color='#2c5fa0', lw=1.4, solid_capstyle='butt')

    # Left branch: dashed → Admissibility
    ax.plot([1.9, 1.9], [6.40, 5.50], color='#808080', lw=1.4,
            linestyle='dashed', solid_capstyle='butt')
    draw_arrow(ax, 1.9, 5.50, 1.9, 5.28, color='#808080', linestyle='dashed')

    # Center branch: solid → P_int
    draw_arrow(ax, 5, 6.40, 5, 5.28, color='#2c5fa0')

    # Right branch: solid → Reduced orbit
    draw_arrow(ax, 8.1, 6.40, 8.1, 5.28, color='#2c5fa0')

    # ── Bottom boxes ────────────────────────────────────────────────────────
    # Left: Admissibility (dashed arrow → topology-robust, χ-independent)
    draw_box(
        ax, 1.9, 4.55,
        "Admissibility\ntopology-robust\n(Sec 6)",
        fc='#dff2df', ec='#307830', fontsize=9.2,
    )

    # Center: P_int diagnostic
    draw_box(
        ax, 5, 4.55,
        "$P_{\\rm int}^{\\rm MX}$\n$C_{\\rm top} = -9\\chi$\n(Sec 4.3)",
        fc='#f5e0ec', ec='#a03065', fontsize=9.2,
    )

    # Right: Reduced orbit
    draw_box(
        ax, 8.1, 4.55,
        "Reduced orbit\n$\\dot{q}^2 + \\chi = 0$\n(Sec 7.1)",
        fc='#ece0f5', ec='#6030a0', fontsize=9.2,
    )

    # ── Legend annotation ───────────────────────────────────────────────────
    ax.plot([0.55, 1.05], [2.80, 2.80], color='#2c5fa0', lw=1.4, label='_')
    ax.annotate('', xy=(1.05, 2.80), xytext=(0.55, 2.80),
                arrowprops=dict(arrowstyle='->', color='#2c5fa0', lw=1.4, mutation_scale=10))
    ax.text(1.15, 2.80, '$\\chi$-controlled', ha='left', va='center', fontsize=8.5)

    ax.plot([0.55, 1.05], [2.35, 2.35], color='#808080', lw=1.4,
            linestyle='dashed', solid_capstyle='butt')
    ax.annotate('', xy=(1.05, 2.35), xytext=(0.55, 2.35),
                arrowprops=dict(arrowstyle='->', color='#808080', lw=1.4,
                               linestyle='dashed', mutation_scale=10))
    ax.text(1.15, 2.35, 'topology-robust ($\\chi$-independent)', ha='left', va='center', fontsize=8.5)

    ax.add_patch(mpatches.FancyBboxPatch(
        (0.40, 2.10), 4.60, 1.00,
        boxstyle='square,pad=0.0', facecolor='none', edgecolor='#aaaaaa', lw=0.8,
    ))

    # ── Save ────────────────────────────────────────────────────────────────
    out = FIGURE_DIR / "fig01_three_level_architecture.png"
    fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
