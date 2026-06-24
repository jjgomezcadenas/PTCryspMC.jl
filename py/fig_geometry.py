#!/usr/bin/env python3
"""Figure: the scanner geometry — a transverse cross-section of the water phantom inside the
crystal ring (segmented into phi blocks), and the (phi, z) block/wheel grid unrolled. Driven by
the real geometry/geometry_sphere.json so the proportions are to scale.

    python3 py/fig_geometry.py            -> docs/figures/geometry.png
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, Rectangle

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "docs", "figures", "geometry.png")

PHANTOM = "#9ecae1"   # water
CRYST = "#d9d9d9"     # crystal
HILITE = "#fc8d62"    # one highlighted block


def main():
    g = json.load(open(os.path.join(REPO, "geometry", "geometry_sphere.json")))
    ph_r = g["phantom"]["radius_cm"]
    sc = g["scanner"]
    r_in = sc["r_inner_cm"]; r_out = r_in + sc["wall_thickness_cm"]
    n_phi = sc["n_phi"]; n_z = sc["n_z"]; H = sc["half_length_cm"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.0))

    # --- (a) transverse cross-section -------------------------------------------------
    dphi = 360.0 / n_phi
    for k in range(n_phi):
        ax1.add_patch(Wedge((0, 0), r_out, k * dphi, (k + 1) * dphi, width=r_out - r_in,
                            facecolor=CRYST, edgecolor="white", linewidth=0.6))
    ax1.add_patch(Wedge((0, 0), r_out, 6 * dphi, 7 * dphi, width=r_out - r_in,
                        facecolor=HILITE, edgecolor="white", linewidth=0.6))
    ax1.add_patch(Circle((0, 0), ph_r, facecolor=PHANTOM, edgecolor="#3182bd", linewidth=1.2))

    ax1.annotate("water phantom\n$R = %.0f$ cm" % ph_r, (0, 0), (0, -ph_r - 9),
                 ha="center", va="top", fontsize=9, color="#08519c",
                 arrowprops=dict(arrowstyle="->", color="#08519c"))
    midr = 0.5 * (r_in + r_out); a = np.deg2rad(6.5 * dphi)
    ax1.annotate("one block\n(a $\\varphi$-segment of\nthe crystal ring)",
                 (midr * np.cos(a), midr * np.sin(a)), (r_out + 4, r_out - 6),
                 ha="left", va="center", fontsize=9, color="#a63603",
                 arrowprops=dict(arrowstyle="->", color="#a63603"))
    ax1.annotate("", (r_in, 0), (0, 0), arrowprops=dict(arrowstyle="<->", color="0.35"))
    ax1.text(0.52 * r_in, 1.6, "$r_\\mathrm{in}=%.1f$ cm" % r_in, fontsize=8.5, color="0.25")
    ax1.text(r_out * 0.62, -r_out * 0.78, "wall %.1f cm" % (r_out - r_in), fontsize=8.5, color="0.35")

    lim = r_out + 12
    ax1.set_xlim(-lim, lim); ax1.set_ylim(-lim, lim); ax1.set_aspect("equal")
    ax1.set_xlabel("x [cm]"); ax1.set_ylabel("y [cm]")
    ax1.set_title("(a)  transverse cross-section  ($%d$ $\\varphi$ blocks)" % n_phi, fontsize=11)

    # --- (b) the (phi, z) block / wheel grid, unrolled --------------------------------
    for k in range(n_phi + 1):
        ax2.plot([k, k], [-H, H], color="0.85", linewidth=0.5)
    for j in range(n_z + 1):
        z = -H + 2 * H * j / n_z
        ax2.plot([0, n_phi], [z, z], color="0.85", linewidth=0.5)
    ax2.add_patch(Rectangle((6, -H + 2 * H * 11 / n_z), 1, 2 * H / n_z,
                            facecolor=HILITE, edgecolor="white"))
    ax2.annotate("one block:\none $(\\varphi, z)$ cell", (6.5, -H + 2 * H * 11.5 / n_z),
                 (n_phi * 0.5, H * 0.55), ha="left", fontsize=9, color="#a63603",
                 arrowprops=dict(arrowstyle="->", color="#a63603"))
    ax2.set_xlim(-1, n_phi + 1); ax2.set_ylim(-H - 4, H + 4)
    ax2.set_xticks([0, n_phi / 2, n_phi]); ax2.set_xticklabels(["0", "$\\pi$", "$2\\pi$"])
    ax2.set_xlabel("azimuth  $\\varphi$  ($%d$ blocks)" % n_phi)
    ax2.set_ylabel("axial  $z$  [cm]  ($%d$ wheels)" % n_z)
    ax2.set_title("(b)  the $(\\varphi, z)$ block / wheel grid", fontsize=11)
    ax2.set_aspect("auto")

    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
