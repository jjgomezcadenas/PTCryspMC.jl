#!/usr/bin/env python3
"""Figure: a coincidence and the line of response (LOR). Left, a true coincidence — two 511 keV
photons leave an annihilation back-to-back, reach two near-opposite blocks, and the LOR through the
two first-interaction points passes through the annihilation. Right, a scatter coincidence — one
photon Compton-scatters in the phantom, so the LOR no longer points at the annihilation.

    python3 py/fig_lor.py            -> docs/figures/lor.png
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "docs", "figures", "lor.png")
PHANTOM = "#9ecae1"; CRYST = "#ececec"; HIT = "#fc8d62"; LOR = "#d7301f"; PH = "#3182bd"


def ray_circle(o, d, R):
    """Smallest positive t with |o + t d| = R."""
    b = 2 * (o[0] * d[0] + o[1] * d[1]); c = o[0]**2 + o[1]**2 - R**2
    disc = b * b - 4 * c
    return (-b + np.sqrt(disc)) / 2.0


def ring(ax, r_in, r_out, n_phi, hits=()):
    dphi = 360.0 / n_phi
    for k in range(n_phi):
        col = HIT if k in hits else CRYST
        ax.add_patch(Wedge((0, 0), r_out, k * dphi, (k + 1) * dphi, width=r_out - r_in,
                           facecolor=col, edgecolor="white", linewidth=0.5))


def block_of(p, n_phi):
    return int((np.degrees(np.arctan2(p[1], p[0])) % 360) / (360.0 / n_phi))


def main():
    g = json.load(open(os.path.join(REPO, "geometry", "geometry_sphere.json")))
    ph_r = g["phantom"]["radius_cm"]; sc = g["scanner"]
    r_in = sc["r_inner_cm"]; r_out = r_in + sc["wall_thickness_cm"]; n_phi = sc["n_phi"]

    X = np.array([2.5, -1.5])                          # annihilation point (inside the phantom)
    d = np.array([np.cos(np.deg2rad(58)), np.sin(np.deg2rad(58))])

    fig, (axT, axS) = plt.subplots(1, 2, figsize=(11, 5.4))

    # --- true coincidence ---
    h1 = X + ray_circle(X, d, r_in) * d
    h2 = X + ray_circle(X, -d, r_in) * (-d)
    ring(axT, r_in, r_out, n_phi, hits=(block_of(h1, n_phi), block_of(h2, n_phi)))
    axT.add_patch(Circle((0, 0), ph_r, facecolor=PHANTOM, edgecolor=PH, lw=1.0, zorder=1))
    for h in (h1, h2):
        axT.annotate("", h, X, arrowprops=dict(arrowstyle="->", color="0.3", lw=1.6))
    axT.plot([h1[0], h2[0]], [h1[1], h2[1]], "--", color=LOR, lw=1.6, zorder=3)
    axT.plot(*X, "x", color="k", ms=9, mew=2, zorder=4)
    for h in (h1, h2):
        axT.plot(*h, "o", color=LOR, ms=6, zorder=4)
    axT.text(X[0] + 1, X[1] - 5, "annihilation", fontsize=9)
    axT.text(0.5 * (h1[0] + h2[0]) - 2, 0.5 * (h1[1] + h2[1]) + 8, "LOR", color=LOR, fontsize=10, rotation=58)
    axT.set_title("(a)  true coincidence\nLOR passes through the annihilation", fontsize=11)

    # --- scatter coincidence: photon 2 Comptons in the phantom ---
    P = X + 4.5 * (-d)                                  # Compton point in the phantom
    d2 = np.array([np.cos(np.deg2rad(205)), np.sin(np.deg2rad(205))])  # kinked direction
    h2s = P + ray_circle(P, d2, r_in) * d2
    ring(axS, r_in, r_out, n_phi, hits=(block_of(h1, n_phi), block_of(h2s, n_phi)))
    axS.add_patch(Circle((0, 0), ph_r, facecolor=PHANTOM, edgecolor=PH, lw=1.0, zorder=1))
    axS.annotate("", h1, X, arrowprops=dict(arrowstyle="->", color="0.3", lw=1.6))
    axS.annotate("", P, X, arrowprops=dict(arrowstyle="->", color="0.3", lw=1.6))
    axS.annotate("", h2s, P, arrowprops=dict(arrowstyle="->", color="0.3", lw=1.6))
    # the reconstructed LOR (between the two hits) and its miss of the annihilation
    t = np.array([h1, h2s])
    axS.plot(t[:, 0], t[:, 1], "--", color=LOR, lw=1.6, zorder=3)
    axS.plot(*X, "x", color="k", ms=9, mew=2, zorder=4)
    axS.plot(*P, "o", color="#6a51a3", ms=7, zorder=4)
    for h in (h1, h2s):
        axS.plot(*h, "o", color=LOR, ms=6, zorder=4)
    axS.text(P[0] - 10, P[1] + 2, "Compton\nscatter", color="#6a51a3", fontsize=9)
    axS.text(X[0] + 1, X[1] + 3, "annihilation", fontsize=9)
    axS.set_title("(b)  scatter coincidence\nLOR misses the annihilation", fontsize=11)

    for ax in (axT, axS):
        lim = r_out + 6
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal"); ax.axis("off")

    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
