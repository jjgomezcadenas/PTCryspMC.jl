#!/usr/bin/env python3
"""Figure: the coincidence-time difference dT = |t1 - t2| of the true coincidences (from a 10^7-
event CsI run), and the fraction of coincidences kept as a function of the window tau. The window
tau is read off the tail of this distribution; the chosen tau = 3 ns keeps ~99% of the pairs.

    python3 py/fig_dt.py            -> docs/figures/dt.png
(needs prod/sphere_water_csi/lors_truth.h5; rerun the production chain to regenerate it)
"""
import os

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "docs", "figures", "dt.png")
LORS = os.path.join(REPO, "prod", "sphere_water_csi", "lors_truth.h5")


def main():
    with h5py.File(LORS) as f:
        dT = np.abs(f["t1_ns"][:] - f["t2_ns"][:])
    tau, dtmax = 3.0, 10.0
    kept = 100.0 * np.mean(dT <= tau)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))

    ax1.hist(dT, bins=200, range=(0, dtmax), color="#4C72B0")
    ax1.axvline(tau, color="#d7301f", ls="--", lw=1.5, label="$\\tau = %.0f$ ns" % tau)
    ax1.set_yscale("log")
    ax1.set_xlabel("$\\Delta T = |t_1 - t_2|$  [ns]")
    ax1.set_ylabel("coincidences / bin")
    ax1.set_title("coincidence-time difference", fontsize=11)
    ax1.legend(frameon=False, fontsize=10)

    s = np.sort(dT); frac = 100.0 * np.arange(1, s.size + 1) / s.size
    m = s <= dtmax
    ax2.plot(s[m], frac[m], color="#55A868", lw=2)
    ax2.axvline(tau, color="#d7301f", ls="--", lw=1.5)
    ax2.axhline(kept, color="#d7301f", ls=":", lw=1)
    ax2.text(tau + 0.2, kept - 6, "$\\tau = %.0f$ ns\nkeeps %.1f%%" % (tau, kept), color="#d7301f", fontsize=9)
    ax2.set_xlim(0, dtmax); ax2.set_ylim(0, 100)
    ax2.set_xlabel("window  $\\tau$  [ns]")
    ax2.set_ylabel("coincidences kept  [%]")
    ax2.set_title("acceptance vs window", fontsize=11)

    fig.suptitle("Choosing the coincidence window (CsI, $10^7$ events)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
