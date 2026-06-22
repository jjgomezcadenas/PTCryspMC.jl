#!/usr/bin/env python3
"""Summary plots for the photon stack written by
scripts/propagate_gammas_in_phantom.jl.

Reads output/phantom_stack.csv (one row per interaction:
event_number, step, x_mm, y_mm, z_mm, e_in_keV, e_dep_keV, process) and draws a
2x2 panel:

  1. escape fraction   — events whose photon leaves the phantom vs. is absorbed
                         inside it (photoelectric / below-cut). The escaping
                         photons are the ones that would reach the detector ring.
  2. scatters/event    — number of Compton scatters, binned 0 / 1 / 2 / >2,
                         over all events.
  3. escape energy     — energy of each escaping photon at the phantom surface
                         (the terminating e_in), escape events only.
  4. escape angle      — that photon's exit direction vs. the initial +z, in
                         degrees, escape events only (large angles = backscatter
                         out the entrance cap).

The transport follows a single photon per event (Compton continues the same
gamma; photoelectric absorbs it; no secondaries), so energy/angle describe that
one photon as it leaves the phantom.

Run from the repo root:
    python3 py/plot_phantom_stack.py
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")  # headless: write a file, never open a window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def per_event(df):
    """Collapse the interaction stack to one row per event."""
    g = df.sort_values(["event_number", "step"]).copy()

    # number of Compton scatters per event
    n_scatter = g.assign(_c=g.process == "compton").groupby("event_number")._c.sum()

    # previous interaction position within the same event (for the last segment)
    for c in ("x_mm", "y_mm", "z_mm"):
        g[c + "_prev"] = g.groupby("event_number")[c].shift(1)

    # the terminating row of each event (largest step)
    last = g.loc[g.groupby("event_number")["step"].idxmax()].set_index("event_number")
    escaped = (last["process"] == "escape").to_numpy()

    # exit energy = energy carried into the terminating interaction
    exit_e = last["e_in_keV"].to_numpy()

    # exit direction = the last travelled segment; for an unscattered photon
    # (single row, no previous point) the direction is the initial +z.
    v = np.stack([last["x_mm"] - last["x_mm_prev"],
                  last["y_mm"] - last["y_mm_prev"],
                  last["z_mm"] - last["z_mm_prev"]], axis=1)
    norm = np.linalg.norm(v, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        cosang = np.where(norm > 0, v[:, 2] / norm, 1.0)
    angle = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))
    angle[last["x_mm_prev"].isna().to_numpy()] = 0.0  # unscattered → 0°

    return pd.DataFrame({
        "n_scatter": n_scatter.astype(int).to_numpy(),
        "escaped": escaped,
        "exit_e": exit_e,
        "angle": angle,
    }, index=last.index)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", default=os.path.join(REPO, "studies", "phantom", "phantom_stack.csv"))
    p.add_argument("--out", default=os.path.join(REPO, "studies", "phantom",
                                                 "phantom_stack_summary.png"))
    args = p.parse_args()

    df = pd.read_csv(args.csv)
    ev = per_event(df)
    n = len(ev)
    esc = ev[ev["escaped"]]                       # escaping photons only
    f_esc = len(esc) / n

    fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
    fig.suptitle(f"Photon stack summary — {n} events  ({os.path.basename(args.csv)})")

    # 1. escape fraction: escaped vs. absorbed inside the phantom
    n_esc = int(ev["escaped"].sum())
    n_abs = n - n_esc
    bars = ax[0, 0].bar(["escaped", "absorbed"], [n_esc, n_abs], color=["C0", "C3"])
    for rect, c in zip(bars, (n_esc, n_abs)):
        ax[0, 0].text(rect.get_x() + rect.get_width() / 2, c,
                      f"{c/n:.1%}", ha="center", va="bottom", fontsize=10)
    ax[0, 0].set_ylabel("events")
    ax[0, 0].set_title(f"Escape fraction = {f_esc:.1%}")

    # 2. scatters per event, binned 0 / 1 / 2 / >2 (all events)
    ns = ev["n_scatter"]
    buckets = [(ns == 0).sum(), (ns == 1).sum(), (ns == 2).sum(), (ns > 2).sum()]
    ax[0, 1].bar(["0", "1", "2", ">2"], buckets, color="C1")
    for i, b in enumerate(buckets):
        ax[0, 1].text(i, b, f"{b/n:.0%}", ha="center", va="bottom", fontsize=9)
    ax[0, 1].set_xlabel("Compton scatters per event")
    ax[0, 1].set_ylabel("events")
    ax[0, 1].set_title("Scatters per event")

    # 3. exit energy of escaping photons
    ax[1, 0].hist(esc["exit_e"], bins=80, range=(0, 520), color="C2")
    ax[1, 0].axvline(511.0, color="gray", ls="--", lw=1.0)
    ax[1, 0].set_xlabel("exit energy [keV]")
    ax[1, 0].set_ylabel("escaping photons")
    ax[1, 0].set_title("Energy of escaping photons")

    # 4. exit angle of escaping photons w.r.t. the initial direction
    ax[1, 1].hist(esc["angle"], bins=90, range=(0, 180), color="C4")
    ax[1, 1].set_xlabel("exit angle to initial +z [deg]")
    ax[1, 1].set_ylabel("escaping photons")
    ax[1, 1].set_title("Angle of escaping photons")

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")
    print(f"  escape fraction = {f_esc:.3f}  ({n_esc} escaped, {n_abs} absorbed)")
    print(f"  scatters 0/1/2/>2 = {buckets}")
    print(f"  unscattered escaping (E>510.9) = {(esc['exit_e'] > 510.9).mean():.3f}")


if __name__ == "__main__":
    main()
