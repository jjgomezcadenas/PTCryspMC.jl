#!/usr/bin/env python3
"""Plot the macroscopic photon cross sections dumped by
scripts/material_xsections.jl.

Reads a CSV with columns energy_keV, compton, phot, pair (sigma in cm^-1) and
writes a log-log plot of each channel plus the total. Pass --material to set the
input/output names and the title in one go, or override with --csv/--out/--title.

Run from the repo root:
    python3 py/plot_material_xsections.py --material Water
    python3 py/plot_material_xsections.py --material CsI
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")  # headless: write a file, never open a window
import matplotlib.pyplot as plt
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--material", default="Water",
                   help="material name (sets default csv/out/title)")
    p.add_argument("--csv", default=None, help="cross-section CSV (default by --material)")
    p.add_argument("--out", default=None, help="output figure (default by --material)")
    p.add_argument("--title", default=None, help="plot title (default by --material)")
    args = p.parse_args()

    slug = args.material.lower()
    csv = args.csv or os.path.join(REPO, "output", f"{slug}_xsections.csv")
    out = args.out or os.path.join(REPO, "output", "control_plots", f"{slug}_xsections.png")
    title = args.title or f"{args.material} photon cross sections"

    df = pd.read_csv(csv)
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
    ax.set_title(title)
    ax.grid(True, which="both", ls=":", alpha=0.4)
    ax.legend()
    fig.tight_layout()

    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
