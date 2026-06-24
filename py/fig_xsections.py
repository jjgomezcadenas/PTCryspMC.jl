#!/usr/bin/env python3
"""Figure: the macroscopic cross sections that drive the transport, Sigma = (mu/rho) * density
[cm^-1], for water (the phantom) and CsI (the crystal), from the NIST XCOM tables in data/. Shows
why 511 keV photons Compton-scatter in the phantom (the patient-scatter background) and why the
crystal both Comptons and photo-absorbs them.

    python3 py/fig_xsections.py            -> docs/figures/xsections.png
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "docs", "figures", "xsections.png")
DENSITY = {"water": 1.0, "CSI": 4.51}    # g/cm^3 (from data/materials.json)


def read_xcom(name):
    """Return (E[MeV], Compton, photoelectric) mass-attenuation columns [cm^2/g]."""
    rows = []
    with open(os.path.join(REPO, "data", f"xcom_{name}.csv")) as f:
        for line in f:
            t = line.split()
            if len(t) == 8 and t[0][0].isdigit():
                rows.append([float(x) for x in t])
    a = np.array(rows)
    return a[:, 0], a[:, 2], a[:, 3]      # energy, incoherent (Compton), photoelectric


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    for ax, (name, title) in zip(axes, [("water", "water phantom"), ("CSI", "CsI crystal")]):
        E, comp, photo = read_xcom(name)
        rho = DENSITY[name]
        ax.loglog(E, comp * rho, color="#1b7837", lw=2, label="Compton")
        ax.loglog(E, photo * rho, color="#762a83", lw=2, label="photoelectric")
        ax.loglog(E, (comp + photo) * rho, color="0.2", lw=2, ls="--", label="total")
        ax.axvline(0.511, color="#d7301f", lw=1.2)
        ax.text(0.55, ax.get_ylim()[0] * 2.0, "511 keV", color="#d7301f", fontsize=9, rotation=90, va="bottom")
        ax.set_xlim(1e-2, 1e1); ax.set_xlabel("photon energy [MeV]")
        ax.set_title(title, fontsize=11); ax.grid(True, which="major", alpha=0.25)
        ax.legend(frameon=False, fontsize=9)
    axes[0].set_ylabel("macroscopic cross section  $\\Sigma$  [cm$^{-1}$]")
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
