#!/usr/bin/env python3
"""Figure: the source activity over the acquisition window. A single-isotope source decays as
A(t) = A_0 e^{-lambda (t - t0)} (15-O, T_half = 2.04 min) over [t0, t1] = [0, 600] s, so the
annihilation times are front-loaded — the origin of the randoms' front-loading.

    python3 py/fig_activity.py            -> docs/figures/activity.png
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "docs", "figures", "activity.png")


def main():
    T_half = 122.24; t0, t1 = 0.0, 600.0
    lam = np.log(2) / T_half
    t = np.linspace(t0, t1, 500)
    A = np.exp(-lam * (t - t0))                         # normalised to A(t0) = 1

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, logy in zip(axes, (False, True)):
        ax.plot(t, A, color="#08519c", lw=2)
        ax.axvline(t0 + T_half, color="0.4", ls="--", lw=1)
        ax.text(t0 + T_half + 10, 0.55, "$T_{1/2}$ = %.0f s" % T_half, color="0.3", fontsize=9)
        ax.set_xlabel("time since acquisition start  $t$  [s]")
        ax.set_ylabel("relative activity  $A(t)/A_0$")
        ax.set_xlim(t0, t1)
        if logy:
            ax.set_yscale("log"); ax.set_title("log scale (a straight line = exponential)", fontsize=11)
        else:
            ax.set_title("linear scale (front-loaded)", fontsize=11)
    fig.suptitle("$^{15}$O activity over the $[0, 600]\\,$s acquisition window", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
