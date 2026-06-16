#!/usr/bin/env python3
"""Diagnostics for the list-mode coincidence file written by
scripts/build_coincidences.jl (truth or detector mode).

Reads a coincidence CSV (one row per coincidence = one LOR):
  event, x1_mm,y1_mm,z1_mm, e1_keV, t1_ns, iz1,iphi1,
         x2_mm,y2_mm,z2_mm, e2_keV, t2_ns, iz2,iphi2,
         x0_mm,y0_mm,z0_mm, truth
and draws a 3x3 panel. This is diagnostics only -- no reconstruction.

  1. energy spectra e1,e2 (+window)   4. ring hit map (phi vs z)   7. source fill x0-z0 (axial)
  2. e1 vs e2                          5. DOI by truth              8. LOR->source distance by truth
  3. e1 by truth (true/scatter)        6. source fill x0-y0 (xverse) 9. truth composition + summary

The LOR->source distance (panel 8) is the perpendicular distance from the true emission
point (x0,y0,z0) to the LOR line through the two hits -- the blur budget (acollinearity +
position resolution + scatter), small for trues, large for scatter.

Run from the repo root (config-driven, recommended):
    python3 py/plot_coincidences.py --config runs/sphere_water_csi.toml
or with explicit paths:
    python3 py/plot_coincidences.py --coinc output/sphere_water_csi/coincidences_det.csv \
                                    --out output/sphere_water_csi/coincidences_det.png
"""
import argparse
import os
import tomllib

import matplotlib

matplotlib.use("Agg")  # headless: write a file, never open a window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def from_config(path):
    """Derive (coinc, out, emin, emax) from a run config TOML — the per-run output dir,
    the detector/truth mode, and the energy-window lines."""
    with open(path, "rb") as f:
        cfg = tomllib.load(f)
    out = cfg.get("output", {})
    tag = out.get("tag") or os.path.splitext(os.path.basename(path))[0]
    rundir = os.path.join(REPO, out.get("dir", "output"), tag)
    det = cfg.get("detector", {})
    eres = det.get("eres", 0.0)
    win = det.get("window_fwhm", 0.0)
    mode = "det" if (det.get("sigma_xyz_mm", 0.0) or eres or win) else "truth"
    coinc = os.path.join(rundir, f"coincidences_{mode}.csv")
    png = os.path.join(rundir, f"coincidences_{mode}.png")
    emin = emax = None
    if eres and win:
        half = win * eres * 511.0
        emin, emax = 511.0 - half, 511.0 + half
    return coinc, png, emin, emax


def lor_source_distance(df):
    """Perpendicular distance [mm] from (x0,y0,z0) to the LOR line through the two hits."""
    P0 = df[["x0_mm", "y0_mm", "z0_mm"]].to_numpy()
    P1 = df[["x1_mm", "y1_mm", "z1_mm"]].to_numpy()
    P2 = df[["x2_mm", "y2_mm", "z2_mm"]].to_numpy()
    u = P2 - P1
    cross = np.cross(P0 - P1, u)
    return np.linalg.norm(cross, axis=1) / np.linalg.norm(u, axis=1)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=None, help="run config TOML (derives coinc/out/window)")
    p.add_argument("--coinc", default=None)
    p.add_argument("--out", default=None)
    # CRYSP1M ring geometry (for DOI and the unrolled hit map)
    p.add_argument("--r-inner-mm", type=float, default=387.0)
    p.add_argument("--wall-mm", type=float, default=37.0)
    p.add_argument("--half-length-mm", type=float, default=512.0)
    # energy-window lines (keV); auto-filled from the config in --config mode
    p.add_argument("--emin", type=float, default=None)
    p.add_argument("--emax", type=float, default=None)
    args = p.parse_args()

    if args.config:
        coinc, out, emin, emax = from_config(args.config)
        args.coinc = args.coinc or coinc
        args.out = args.out or out
        if args.emin is None: args.emin = emin
        if args.emax is None: args.emax = emax
    if not args.coinc:
        p.error("pass --config or --coinc")
    args.out = args.out or os.path.splitext(args.coinc)[0] + ".png"

    Ri, wall, H = args.r_inner_mm, args.wall_mm, args.half_length_mm
    df = pd.read_csv(args.coinc)
    n = len(df)
    ist = (df["truth"] == "true").to_numpy()
    isc = ~ist

    # both-hit arrays (panels that use every detected hit), tagged with the pair's truth
    e_all = np.concatenate([df["e1_keV"], df["e2_keV"]])
    x_all = np.concatenate([df["x1_mm"], df["x2_mm"]])
    y_all = np.concatenate([df["y1_mm"], df["y2_mm"]])
    z_all = np.concatenate([df["z1_mm"], df["z2_mm"]])
    true_all = np.concatenate([ist, ist])
    r_all = np.hypot(x_all, y_all)
    doi_all = r_all - Ri
    phi_all = np.degrees(np.arctan2(y_all, x_all) % (2 * np.pi))

    d_lor = lor_source_distance(df)

    fig = plt.figure(figsize=(15, 13))
    fig.suptitle(f"Coincidences — {n} LORs ({os.path.basename(args.coinc)})", fontsize=13)

    def ax(i):
        return fig.add_subplot(3, 3, i)

    # 1. energy spectra e1, e2 with the window + 511 line
    a = ax(1)
    a.hist(df["e1_keV"], bins=80, range=(0, 560), histtype="step", color="C0", label="e1")
    a.hist(df["e2_keV"], bins=80, range=(0, 560), histtype="step", color="C1", label="e2")
    a.axvline(511, color="gray", ls="--", lw=1)
    for ev in (args.emin, args.emax):
        if ev is not None:
            a.axvline(ev, color="C3", ls=":", lw=1)
    a.set_yscale("log")
    a.set_xlabel("energy [keV]"); a.set_ylabel("hits"); a.legend()
    a.set_title("Energy spectra (e1, e2)")

    # 2. e1 vs e2
    a = ax(2)
    h = a.hist2d(df["e1_keV"], df["e2_keV"], bins=70, range=[[250, 560], [250, 560]], cmap="viridis")
    fig.colorbar(h[3], ax=a, label="LORs")
    a.set_xlabel("e1 [keV]"); a.set_ylabel("e2 [keV]"); a.set_aspect("equal")
    a.set_title("e1 vs e2 (photopeak at 511,511)")

    # 3. e1 split by truth
    a = ax(3)
    a.hist(df.loc[ist, "e1_keV"], bins=80, range=(0, 560), histtype="step", color="C2", label="true")
    a.hist(df.loc[isc, "e1_keV"], bins=80, range=(0, 560), histtype="step", color="C1", label="scatter")
    a.axvline(511, color="gray", ls="--", lw=1)
    a.set_yscale("log")
    a.set_xlabel("e1 [keV]"); a.set_ylabel("LORs"); a.legend()
    a.set_title("Energy by truth")

    # 4. ring hit map (unrolled)
    a = ax(4)
    h = a.hist2d(phi_all, z_all, bins=[72, 60], range=[[0, 360], [-H, H]], cmap="viridis")
    fig.colorbar(h[3], ax=a, label="hits")
    a.set_xlabel("φ [deg]"); a.set_ylabel("z [mm]")
    a.set_title("Ring hit map (unrolled)")

    # 5. DOI by truth
    a = ax(5)
    a.hist(doi_all[true_all], bins=50, range=(0, wall), histtype="step", color="C2", label="true")
    a.hist(doi_all[~true_all], bins=50, range=(0, wall), histtype="step", color="C1", label="scatter")
    a.set_xlabel("DOI = r − R_inner [mm]"); a.set_ylabel("hits"); a.legend()
    a.set_title("Depth of interaction")

    # 6. source fill, transverse
    a = ax(6)
    lim = 1.05 * max(np.abs(df["x0_mm"]).max(), np.abs(df["y0_mm"]).max())
    h = a.hist2d(df["x0_mm"], df["y0_mm"], bins=60, range=[[-lim, lim], [-lim, lim]], cmap="magma")
    fig.colorbar(h[3], ax=a, label="emissions")
    a.set_xlabel("x0 [mm]"); a.set_ylabel("y0 [mm]"); a.set_aspect("equal")
    a.set_title("Source fill (transverse)")

    # 7. source fill, axial
    a = ax(7)
    limz = 1.05 * np.abs(df["z0_mm"]).max()
    h = a.hist2d(df["x0_mm"], df["z0_mm"], bins=60, range=[[-lim, lim], [-limz, limz]], cmap="magma")
    fig.colorbar(h[3], ax=a, label="emissions")
    a.set_xlabel("x0 [mm]"); a.set_ylabel("z0 [mm]")
    a.set_title("Source fill (axial)")

    # 8. LOR -> source distance by truth (the blur budget)
    a = ax(8)
    hi = np.percentile(d_lor, 99)
    a.hist(d_lor[ist], bins=60, range=(0, hi), histtype="step", color="C2", label="true")
    a.hist(d_lor[isc], bins=60, range=(0, hi), histtype="step", color="C1", label="scatter")
    a.set_yscale("log")
    a.set_xlabel("LOR → emission-point distance [mm]"); a.set_ylabel("LORs"); a.legend()
    a.set_title("LOR localization (blur budget)")

    # 9. truth composition + summary text
    a = ax(9)
    n_true, n_sc = int(ist.sum()), int(isc.sum())
    bars = a.bar(["true", "scatter"], [n_true, n_sc], color=["C2", "C1"])
    for b, c in zip(bars, (n_true, n_sc)):
        a.text(b.get_x() + b.get_width() / 2, c, f"{c/max(n,1):.1%}", ha="center", va="bottom")
    a.set_ylabel("LORs"); a.set_title("Truth composition")
    txt = (f"n LORs = {n}\n"
           f"true = {n_true} ({n_true/max(n,1):.1%})\n"
           f"scatter = {n_sc} ({n_sc/max(n,1):.1%})\n"
           f"median LOR→src: true {np.median(d_lor[ist]):.1f} mm\n"
           f"                scatter {np.median(d_lor[isc]):.1f} mm")
    a.text(0.97, 0.97, txt, transform=a.transAxes, ha="right", va="top", fontsize=9,
           bbox=dict(boxstyle="round", fc="white", ec="0.7"))

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=120)
    print(f"wrote {args.out}")
    print(f"  {n} LORs;  true {n_true/max(n,1):.1%} / scatter {n_sc/max(n,1):.1%};  "
          f"median LOR→src (true) {np.median(d_lor[ist]):.1f} mm")


if __name__ == "__main__":
    main()
