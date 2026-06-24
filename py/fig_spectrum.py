#!/usr/bin/env python3
"""Figure: the single-photon energy spectrum measured in the crystal (CsI, 10^7 events), smeared by
the 5% energy resolution. The full-absorption photopeak sits at 511 keV; below it is the Compton
continuum (photons that deposited only the recoil and escaped), with its edge at ~341 keV. The two
lower energy cuts are marked: the studies cut (300 keV, keeps the Compton shoulder) and the
reconstruction cut (450 keV, selects the photopeak).

    python3 py/fig_spectrum.py            -> docs/figures/spectrum.png
(needs prod/sphere_water_csi/singles.h5)
"""
import os

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "docs", "figures", "spectrum.png")
SINGLES = os.path.join(REPO, "prod", "sphere_water_csi", "singles.h5")
FWHM2SIG = 2.3548200450309493


def main():
    with h5py.File(SINGLES) as f:
        e = f["e_keV"][:].astype(np.float64) * 0.1     # decode Int16 @ 0.1 keV
    e = e[e > 0]
    rng = np.random.default_rng(0)                     # 5% FWHM @ 511 keV detector smearing
    sigma = 0.05 * np.sqrt(511.0 * e) / FWHM2SIG
    es = e + sigma * rng.standard_normal(e.size)

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.hist(es, bins=325, range=(0, 650), color="#bdbdbd", edgecolor="none")
    ax.axvline(511, color="#08519c", lw=1.2)
    ax.text(513, ax.get_ylim()[1] * 0.92, "511 keV\nphotopeak", color="#08519c", fontsize=9)
    ax.axvline(340.7, color="#1b7837", lw=1.0, ls=":")
    ax.text(300, ax.get_ylim()[1] * 0.55, "Compton\nedge", color="#1b7837", fontsize=9, ha="right")
    for x, lab, col in [(300, "studies cut\n300 keV", "#762a83"), (450, "reco cut\n450 keV", "#d7301f")]:
        ax.axvline(x, color=col, lw=1.6, ls="--")
        ax.text(x - 5, ax.get_ylim()[1] * 0.25, lab, color=col, fontsize=9, ha="right", rotation=90, va="center")
    ax.set_xlim(0, 650); ax.set_xlabel("deposited energy (smeared)  [keV]")
    ax.set_ylabel("singles / bin")
    ax.set_title("Crystal energy spectrum and the energy cuts (CsI)", fontsize=11)
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
