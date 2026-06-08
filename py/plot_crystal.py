#!/usr/bin/env python3
"""Containment / topology study for a scintillator crystal, from the stack
written by scripts/shoot_gammas_to_crystal.jl (CsI, BGO, ...).

Reads a crystal stack CSV (one row per interaction:
event_number, step, x_mm, y_mm, z_mm, e_in_keV, e_dep_keV, process) and draws a
3x3 panel:

  1. containment       — photon fully absorbed in the crystal vs. escaped
  2. photo (contained) — contained events with 0 Compton (direct photoelectric)
  3. compton (contained) — contained events with >=1 Compton
  4. Etot (contained)  — energy deposited per contained event (the 511 peak)
  5. scatters (contained) — Compton scatters per contained event, 0/1/2/>2
  6. z1 (contained)    — depth of the 1st Compton interaction
  7. z2 (contained)    — depth of the 2nd Compton interaction
  8. |z1 - z2| (contained) — first-to-second Compton depth separation
  9. x vs y (contained) — transverse map of all interaction points

"contained" = the photon did not escape (terminating process != escape).
All panels but #1 are restricted to contained events. z is the depth from the
entry face (the crystal spans z in [0, depth]).

Run from the repo root:
    python3 py/plot_crystal.py --csv output/bgo_crystal_cone_stack.csv \
                               --out output/control_plots/bgo_crystal_cone_summary.png
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")  # headless: write a file, never open a window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", default=os.path.join(REPO, "output", "csi_crystal_stack.csv"))
    p.add_argument("--out", default=os.path.join(REPO, "output", "control_plots",
                                                 "csi_crystal_summary.png"))
    args = p.parse_args()

    df = pd.read_csv(args.csv).sort_values(["event_number", "step"])
    grp = df.groupby("event_number")

    # per-event quantities
    contained = (grp["process"].last() != "escape")          # bool, indexed by event
    n_scatter = (df["process"] == "compton").groupby(df["event_number"]).sum()
    etot = grp["e_dep_keV"].sum()
    ev_contained = contained.index[contained]                # contained event ids
    n = len(contained)
    n_cont = int(contained.sum())

    # photo / compton among contained
    nsc_c = n_scatter.loc[ev_contained]
    n_photo = int((nsc_c == 0).sum())
    n_compton = int((nsc_c >= 1).sum())

    # Compton interactions in contained events, ranked within each event
    comp = df[(df["process"] == "compton") & df["event_number"].isin(ev_contained)].copy()
    comp["rank"] = comp.groupby("event_number").cumcount()
    z1 = comp.loc[comp["rank"] == 0].set_index("event_number")["z_mm"]
    z2 = comp.loc[comp["rank"] == 1].set_index("event_number")["z_mm"]
    dz = (z1 - z2).abs().dropna()

    # all interaction points in contained events (no escape rows there)
    inter = df[df["event_number"].isin(ev_contained)]

    depth = df["z_mm"].max()
    half = df[["x_mm", "y_mm"]].abs().max().max()

    fig, ax = plt.subplots(3, 3, figsize=(14, 12))
    fig.suptitle(f"Crystal study — {n} events  ({os.path.basename(args.csv)})", fontsize=13)

    def frac_bar(a, labels, counts, denom, title):
        bars = a.bar(labels, counts, color=["C0", "C7"])
        for r, c in zip(bars, counts):
            a.text(r.get_x() + r.get_width() / 2, c, f"{c/denom:.1%}",
                   ha="center", va="bottom", fontsize=10)
        a.set_ylabel("events")
        a.set_title(title)

    # 1. containment
    frac_bar(ax[0, 0], ["contained", "escaped"], [n_cont, n - n_cont], n,
             f"Containment = {n_cont/n:.1%}")
    # 2. photo (of contained)
    frac_bar(ax[0, 1], ["photo", "other"], [n_photo, n_cont - n_photo], n_cont,
             f"Photo (of contained) = {n_photo/n_cont:.1%}")
    # 3. compton (of contained)
    frac_bar(ax[0, 2], ["compton", "other"], [n_compton, n_cont - n_compton], n_cont,
             f"Compton (of contained) = {n_compton/n_cont:.1%}")

    # 4. Etot (contained)
    ax[1, 0].hist(etot.loc[ev_contained], bins=80, range=(0, 520), color="C2")
    ax[1, 0].axvline(511.0, color="gray", ls="--", lw=1.0)
    ax[1, 0].set_xlabel("Etot [keV]"); ax[1, 0].set_ylabel("contained events")
    ax[1, 0].set_title("Energy deposited (contained)")

    # 5. scatters per contained event
    buckets = [(nsc_c == 0).sum(), (nsc_c == 1).sum(), (nsc_c == 2).sum(), (nsc_c > 2).sum()]
    bb = ax[1, 1].bar(["0", "1", "2", ">2"], buckets, color="C1")
    for r, c in zip(bb, buckets):
        ax[1, 1].text(r.get_x() + r.get_width() / 2, c, f"{c/n_cont:.0%}",
                      ha="center", va="bottom", fontsize=9)
    ax[1, 1].set_xlabel("Compton scatters"); ax[1, 1].set_ylabel("contained events")
    ax[1, 1].set_title("Scatters per event (contained)")

    # 6. z1
    ax[1, 2].hist(z1, bins=60, range=(0, depth), color="C3")
    ax[1, 2].set_xlabel("z1 [mm]"); ax[1, 2].set_ylabel("events")
    ax[1, 2].set_title("Depth of 1st Compton (contained)")

    # 7. z2
    ax[2, 0].hist(z2, bins=60, range=(0, depth), color="C4")
    ax[2, 0].set_xlabel("z2 [mm]"); ax[2, 0].set_ylabel("events")
    ax[2, 0].set_title("Depth of 2nd Compton (contained)")

    # 8. |z1 - z2|
    ax[2, 1].hist(dz, bins=60, range=(0, depth), color="C5")
    ax[2, 1].set_xlabel("|z1 - z2| [mm]"); ax[2, 1].set_ylabel("events")
    ax[2, 1].set_title("1st-to-2nd Compton separation")

    # 9. x vs y (all interactions in contained events)
    h = ax[2, 2].hist2d(inter["x_mm"], inter["y_mm"], bins=60,
                        range=[[-half, half], [-half, half]], cmap="viridis")
    fig.colorbar(h[3], ax=ax[2, 2], label="interactions")
    ax[2, 2].set_xlabel("x [mm]"); ax[2, 2].set_ylabel("y [mm]")
    ax[2, 2].set_title("Interaction map (contained)")
    ax[2, 2].set_aspect("equal")

    fig.tight_layout(rect=(0, 0, 1, 0.98))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=130)
    print(f"wrote {args.out}")
    print(f"  contained = {n_cont/n:.3f}  (photo {n_photo/n_cont:.3f}, "
          f"compton {n_compton/n_cont:.3f} of contained)")
    print(f"  events with >=2 Compton (z2 defined) = {len(dz)}")


if __name__ == "__main__":
    main()
