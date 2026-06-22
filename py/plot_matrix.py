#!/usr/bin/env python3
"""Cross-config comparison of the list-mode coincidence files — one figure putting the
run conditions side by side. Diagnostics only, no reconstruction.

Reads N run configs (default: all runs/*.toml), each pointing at its
output/<tag>/coincidences_{truth,det}.csv, and draws a 2x2 comparison:

  1. energy spectra overlaid (normalised)   2. true / scatter fraction per config
  3. coincidence acceptance per config       4. LOR resolution (median LOR->source, true)

Run from the repo root (after scripts/run_matrix.sh has produced the coincidence files):
    python3 py/plot_matrix.py                                  # all runs/*.toml
    python3 py/plot_matrix.py sphere_water_csi sphere_water_bgo
"""
import argparse
import glob
import os
import sys
import tomllib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_coincidences import from_config, lor_source_distance  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def resolve(c):
    """A config arg (name, name.toml, or a path) -> runs/<name>.toml under the repo."""
    if os.path.isabs(c) or os.path.exists(c):
        return c
    return os.path.join(REPO, "runs", os.path.splitext(os.path.basename(c))[0] + ".toml")


def load_run(cfgpath):
    coinc, _, _, _ = from_config(cfgpath)
    with open(cfgpath, "rb") as f:
        cfg = tomllib.load(f)
    nev = cfg.get("source", {}).get("nevents", None)
    tag = os.path.splitext(os.path.basename(cfgpath))[0]
    df = pd.read_csv(coinc)
    df["truth_is"] = df["truth"] == "true"
    df["d_lor"] = lor_source_distance(df)
    return {"tag": tag, "nev": nev, "df": df}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("configs", nargs="*", help="run configs (default: all runs/*.toml)")
    p.add_argument("--out", default=os.path.join(REPO, "output", "control_plots", "matrix_summary.png"))
    args = p.parse_args()

    cfgs = [resolve(c) for c in args.configs] if args.configs \
        else sorted(glob.glob(os.path.join(REPO, "runs", "*.toml")))

    runs = []
    for c in cfgs:
        try:
            runs.append(load_run(c))
        except FileNotFoundError:
            print(f"skip {os.path.basename(c)}: no coincidence file (run build first)")
    if not runs:
        sys.exit("no coincidence files found — run scripts/run_matrix.sh first")

    tags = [r["tag"] for r in runs]
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(runs), 1)))
    y = np.arange(len(runs))

    fig, axs = plt.subplots(2, 2, figsize=(15, 11))
    fig.suptitle(f"Coincidence matrix — {len(runs)} configs", fontsize=14)

    # 1. energy spectra overlaid (both hits, normalised so shapes compare)
    a = axs[0, 0]
    for r, c in zip(runs, colors):
        e = np.concatenate([r["df"]["e1_keV"], r["df"]["e2_keV"]])
        a.hist(e, bins=80, range=(250, 560), histtype="step", density=True,
               color=c, label=r["tag"])
    a.axvline(511, color="gray", ls="--", lw=1)
    a.set_yscale("log")
    a.set_xlabel("energy [keV]"); a.set_ylabel("normalised hits")
    a.set_title("Energy spectra (Compton shoulder → photopeak)")
    a.legend(fontsize=8)

    # 2. true / scatter fraction per config (stacked horizontal bars)
    a = axs[0, 1]
    tf = np.array([r["df"]["truth_is"].mean() for r in runs])
    a.barh(y, tf, color="C2", label="true")
    a.barh(y, 1 - tf, left=tf, color="C1", label="scatter")
    for i, f in enumerate(tf):
        a.text(0.02, i, f"{f:.0%}", va="center", fontsize=8, color="white")
    a.set_yticks(y); a.set_yticklabels(tags, fontsize=8); a.invert_yaxis()
    a.set_xlabel("fraction"); a.set_xlim(0, 1); a.legend(loc="lower right", fontsize=8)
    a.set_title("True / scatter")

    # 3. coincidence acceptance per config
    a = axs[1, 0]
    acc = [100 * len(r["df"]) / r["nev"] if r["nev"] else np.nan for r in runs]
    bars = a.barh(y, acc, color=colors)
    for i, (v, r) in enumerate(zip(acc, runs)):
        a.text(v, i, f" {v:.1f}% ({len(r['df'])})", va="center", fontsize=8)
    a.set_yticks(y); a.set_yticklabels(tags, fontsize=8); a.invert_yaxis()
    a.set_xlabel("coincidence acceptance [% of annihilations]")
    a.set_title("Acceptance (yield)")

    # 4. LOR resolution: median LOR->source for trues
    a = axs[1, 1]
    med = [np.median(r["df"].loc[r["df"]["truth_is"], "d_lor"]) if r["df"]["truth_is"].any()
           else np.nan for r in runs]
    a.barh(y, med, color=colors)
    for i, v in enumerate(med):
        a.text(v, i, f" {v:.1f} mm", va="center", fontsize=8)
    a.set_yticks(y); a.set_yticklabels(tags, fontsize=8); a.invert_yaxis()
    a.set_xlabel("median LOR → emission-point distance, true [mm]")
    a.set_title("LOR resolution (blur budget)")

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=120)
    print(f"wrote {args.out}  ({len(runs)} configs: {', '.join(tags)})")


if __name__ == "__main__":
    main()
