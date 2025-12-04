#!/usr/bin/env python3
"""
Generate a simple schematic for the response chain:
Pi -> C_Pi -> tau_G -> chi(omega) -> G_eff(Box) -> R_{mu nu}

Outputs:
- /home/ovidiu/EESH_stablized/TRIAD_Geometry/Latex/fig/chain_overview.png
"""
import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch


def add_box(ax, xy, text, boxstyle="round,pad=0.3,rounding_size=0.05"):
    x, y = xy
    width, height = 0.26, 0.16
    rect = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle=boxstyle,
        linewidth=1.2,
        edgecolor="black",
        facecolor="#f5f5f5",
    )
    ax.add_patch(rect)
    ax.text(x, y, text, ha="center", va="center", fontsize=10)
    return rect


def add_arrow(ax, p0, p1):
    arrow = FancyArrowPatch(
        p0,
        p1,
        arrowstyle="-|>",
        mutation_scale=12,
        linewidth=1.2,
        color="black",
        shrinkA=8,
        shrinkB=8,
    )
    ax.add_patch(arrow)


def main():
    out_path = "/home/ovidiu/EESH_stablized/TRIAD_Geometry/Latex/fig/chain_overview.png"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 2.4), constrained_layout=True)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    # X positions for nodes
    xs = [0.07, 0.24, 0.41, 0.58, 0.75, 0.92]
    y = 0.5

    nodes = [
        "Pi = sqrt(F_Q) · S · I",
        "C_Pi(τ)",
        "tau_G = ∫ C_Pi",
        "chi(ω) ≈ 1 - i ω tau_G",
        "G_eff(□) = G0 · chi(□)",
        "R_{μν}",
    ]

    # Draw boxes
    for x, label in zip(xs, nodes):
        add_box(ax, (x, y), label)

    # Draw arrows
    for i in range(len(xs) - 1):
        add_arrow(ax, (xs[i] + 0.06, y), (xs[i + 1] - 0.06, y))

    # Footer (assumptions band)
    ax.text(
        0.5,
        0.12,
        r"Band: $|\omega|\,\tau_G\ll 1$, LWSS+mixing, passivity; DC via calibration",
        ha="center",
        va="center",
        fontsize=9,
        color="#444444",
    )

    fig.savefig(out_path, dpi=200)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()


