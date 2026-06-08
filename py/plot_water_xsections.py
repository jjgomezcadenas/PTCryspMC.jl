#!/usr/bin/env python3
"""Plot the macroscopic photon cross sections dumped by
scripts/water_xsections.jl.

Reads a CSV with columns energy_keV, compton, phot, pair (sigma in cm^-1) and
writes a log-log plot of each channel plus the total to output/control_plots/.

Run from the repo root:
    python3 py/plot_water_xsections.py
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")  # headless: write a file, never open a window
import matplotlib.pyplot as plt
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv", default=os.path.join(REPO, "output", "water_xsections.csv"),
                   help="cross-section CSV from water_xsections.jl")
    p.add_argument("--out", default=os.path.join(REPO, "output", "control_plots",
                                                 "water_xsections.png"),
                   help="output figure")
    p.add_argument("--title", default="Water photon cross sections")
    args = p.parse_args()

    df = pd.read_csv(args.csv)
    total = df["compton"] + df["phot"] + df["pair"]

    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    ax.loglog(df["energy_keV"], total,        "-",  color="black", lw=2.0, label="total")
    ax.loglog(df["energy_keV"], df["compton"], "o-", ms=4, label="Compton")
    ax.loglog(df["energy_keV"], df["phot"],    "s-", ms=4, label="photoelectric")
    # Pair is zero below threshold; log axis can't show zeros, so plot only > 0.
    pair_on = df[df["pair"] > 0.0]
    ax.loglog(pair_on["energy_keV"], pair_on["pair"], "^-", ms=4, label="pair")

    ax.axvline(511.0, color="gray", ls="--", lw=1.0)
    ax.text(511.0, ax.get_ylim()[0], " 511 keV", color="gray",
            va="bottom", ha="left", fontsize=8)

    ax.set_xlabel("photon energy [keV]")
    ax.set_ylabel(r"macroscopic cross section $\Sigma$ [cm$^{-1}$]")
    ax.set_title(args.title)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend()
    fig.tight_layout()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
