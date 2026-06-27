#!/usr/bin/env python3
"""Control plots for a production LOR file (prod/<tag>/lors_det.h5).

A 3x3 panel of QA diagnostics for the list-mode deliverable — energy, the ring hit map,
the LOR axial midpoint, the source fill, the coincidence Δt, and the truth/scatter
composition. Reads ONLY lors_det.h5 (so it runs on already-pruned runs) plus the scanner
geometry from the config's geometry JSON. Writes prod/<tag>/<tag>_lors_det.png.

This is the production counterpart of py/plot_coincidences.py (which is the dev-chain
plotter: CSV from output/, two truth classes). Here the input is HDF5, truth has three
classes (true/scatter/random), and the run provenance comes from the HDF5 attrs.

Run from the repo root (config-driven):
    python3 py/plot_prod.py --config runs/sphere_water_bgo_1MBq.toml
or with explicit paths:
    python3 py/plot_prod.py --lors prod/<tag>/lors_det.h5 --out prod/<tag>/lors_det.png
"""
import argparse
import json
import os
import tomllib

import h5py
import matplotlib

matplotlib.use("Agg")  # headless: write a file, never open a window
import matplotlib.pyplot as plt
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Truth-class colors, consistent across every panel (colorblind-safe Tableau).
C_TRUE, C_SCAT, C_RAND = "#2ca02c", "#ff7f0e", "#d62728"
C_SINGLE, C_MULTI = "#fdae6b", "#e6550d"   # single vs multiple scatter
TRUTH_TRUE, TRUTH_SCATTER, TRUTH_RANDOM = 0, 1, 2


def paths_from_config(path):
    """(lors, out, geom_file) from a run config TOML — the production output dir + the
    geometry JSON for the scanner dimensions."""
    with open(path, "rb") as f:
        cfg = tomllib.load(f)
    out = cfg.get("output", {})
    tag = out.get("tag") or os.path.splitext(os.path.basename(path))[0]
    prod_dir = out.get("prod_dir", "prod")
    rundir = os.path.join(REPO, prod_dir, tag)
    geom_file = os.path.join(REPO, cfg.get("geometry", {}).get("file", "geometry/geometry.json"))
    # The PNG carries the full tag (the h5 keeps its short name): control plots from several
    # runs get opened side by side, so the filename must identify the run on its own.
    return os.path.join(rundir, "lors_det.h5"), os.path.join(rundir, f"{tag}_lors_det.png"), geom_file


def scanner_geometry(geom_file):
    """Ring (r_inner_mm, wall_mm, half_length_mm) from the geometry JSON's scanner section."""
    with open(geom_file) as f:
        g = json.load(f)
    s = g["scanner"]
    return 10.0 * s["r_inner_cm"], 10.0 * s["wall_thickness_cm"], 10.0 * s["half_length_cm"]


def load_lors(lors_file):
    """Read lors_det.h5 into a dict of float arrays (de-quantizing the Int16 columns via the
    stored scales) plus the provenance attrs."""
    with h5py.File(lors_file, "r") as f:
        a = dict(f.attrs)
        xs = float(a.get("xyz_scale_mm", 0.1))
        es = float(a.get("e_scale_keV", 0.1))
        col = lambda k, s=1.0: f[k][:].astype(np.float64) * s
        d = {}
        for k in ("x1_mm", "y1_mm", "z1_mm", "x2_mm", "y2_mm", "z2_mm", "x0_mm", "y0_mm", "z0_mm"):
            d[k] = col(k, xs)
        for k in ("e1_keV", "e2_keV"):
            d[k] = col(k, es)
        for k in ("t1_ns", "t2_ns", "dt_ns"):
            d[k] = col(k)
        for k in ("truth", "nscat1", "nscat2"):
            d[k] = f[k][:].astype(np.int64)
    return d, a


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=None, help="run config TOML (derives lors/out/geometry)")
    p.add_argument("--lors", default=None, help="path to lors_det.h5")
    p.add_argument("--out", default=None, help="output PNG path")
    p.add_argument("--geometry", default=None, help="geometry JSON (scanner dimensions)")
    args = p.parse_args()

    if args.config:
        lors, out, geom = paths_from_config(args.config)
        args.lors = args.lors or lors
        args.out = args.out or out
        args.geometry = args.geometry or geom
    if not args.lors:
        p.error("pass --config or --lors")
    if not args.out:                                  # tag the PNG by its run dir (prod/<tag>/)
        rundir = os.path.dirname(os.path.abspath(args.lors))
        args.out = os.path.join(rundir, f"{os.path.basename(rundir)}_lors_det.png")
    Ri, wall, H = scanner_geometry(args.geometry) if args.geometry else (387.0, 37.0, 512.0)

    d, attrs = load_lors(args.lors)
    truth = d["truth"]
    n = truth.size
    is_t = truth == TRUTH_TRUE
    is_s = truth == TRUTH_SCATTER
    is_r = truth == TRUTH_RANDOM
    n_t, n_s, n_r = int(is_t.sum()), int(is_s.sum()), int(is_r.sum())
    nev = int(attrs.get("nevents", 0))
    crystal = attrs.get("crystal", b"?")
    crystal = crystal.decode() if isinstance(crystal, bytes) else crystal
    tag = attrs.get("scenario_tag", b"?")
    tag = tag.decode() if isinstance(tag, bytes) else tag
    eres = float(attrs.get("eres", 0.0))
    emin = float(attrs.get("emin_keV", 0.0))
    tau = float(attrs.get("tau_ns", 0.0))
    sig = float(attrs.get("sigma_xyz_mm", 0.0))

    # both-hit arrays (per-hit panels), tagged with each pair's truth
    e_all = np.concatenate([d["e1_keV"], d["e2_keV"]])
    x_all = np.concatenate([d["x1_mm"], d["x2_mm"]])
    y_all = np.concatenate([d["y1_mm"], d["y2_mm"]])
    z_all = np.concatenate([d["z1_mm"], d["z2_mm"]])
    truth_all = np.concatenate([truth, truth])
    r_all = np.hypot(x_all, y_all)
    doi_all = r_all - Ri
    phi_all = np.degrees(np.arctan2(y_all, x_all) % (2 * np.pi))
    zmid = 0.5 * (d["z1_mm"] + d["z2_mm"])                 # LOR axial midpoint
    nsc = d["nscat1"] + d["nscat2"]                        # total scatter multiplicity

    fig = plt.figure(figsize=(15, 13))
    fig.suptitle(f"{tag} — {n:,} LORs ({crystal}, eres {eres:.0%}, Emin {emin:.0f} keV, τ {tau:g} ns)",
                 fontsize=13)
    ax = lambda i: fig.add_subplot(3, 3, i)

    # 1. energy spectra e1, e2 + cut + 511
    a = ax(1)
    a.hist(d["e1_keV"], bins=80, range=(0, 560), histtype="step", color="C0", label="e1")
    a.hist(d["e2_keV"], bins=80, range=(0, 560), histtype="step", color="C1", label="e2")
    a.axvline(511, color="gray", ls="--", lw=1)
    if emin:
        a.axvline(emin, color=C_RAND, ls=":", lw=1.2, label=f"Emin {emin:.0f}")
    a.set_yscale("log"); a.set_xlabel("energy [keV]"); a.set_ylabel("hits"); a.legend()
    a.set_title("Energy spectra (e1, e2)")

    # 2. e1 vs e2
    a = ax(2)
    h = a.hist2d(d["e1_keV"], d["e2_keV"], bins=70, range=[[300, 560], [300, 560]], cmap="viridis")
    fig.colorbar(h[3], ax=a, label="LORs")
    a.set_xlabel("e1 [keV]"); a.set_ylabel("e2 [keV]"); a.set_aspect("equal")
    a.set_title("e1 vs e2 (photopeak 511,511)")

    # 3. energy by truth class
    a = ax(3)
    for mask, c, lab in ((truth_all == TRUTH_TRUE, C_TRUE, "true"),
                         (truth_all == TRUTH_SCATTER, C_SCAT, "scatter"),
                         (truth_all == TRUTH_RANDOM, C_RAND, "random")):
        if mask.any():
            a.hist(e_all[mask], bins=80, range=(0, 560), histtype="step", color=c, label=lab)
    a.axvline(511, color="gray", ls="--", lw=1)
    a.set_yscale("log"); a.set_xlabel("energy [keV]"); a.set_ylabel("hits"); a.legend()
    a.set_title("Energy by truth")

    # 4. ring hit map (unrolled φ–z)
    a = ax(4)
    h = a.hist2d(phi_all, z_all, bins=[72, 40], range=[[0, 360], [-H, H]], cmap="viridis")
    fig.colorbar(h[3], ax=a, label="hits")
    a.set_xlabel("φ [deg]"); a.set_ylabel("z [mm]")
    a.set_title("Ring hit map (unrolled)")

    # 5. axial efficiency / sensitivity profile (LOR midpoint z)
    a = ax(5)
    a.hist(zmid, bins=60, range=(-H, H), histtype="step", color="C0", label="all")
    if n_t:
        a.hist(zmid[is_t], bins=60, range=(-H, H), histtype="step", color=C_TRUE, label="true")
    a.set_xlabel("LOR axial midpoint z [mm]"); a.set_ylabel("LORs"); a.legend()
    # NB: the midpoint of a back-to-back pair tracks the SOURCE, not the ring — for a compact
    # central source this traces the source axial extent, not the detector sensitivity (which would
    # need a FOV-filling source). The full ring length shows up in the hit map (panel 4) instead.
    a.set_title("LOR axial midpoint (≈ source axial extent)")

    # 6. source fill, transverse
    a = ax(6)
    lim = 1.05 * max(np.abs(d["x0_mm"]).max(), np.abs(d["y0_mm"]).max(), 1.0)
    h = a.hist2d(d["x0_mm"], d["y0_mm"], bins=60, range=[[-lim, lim], [-lim, lim]], cmap="magma")
    fig.colorbar(h[3], ax=a, label="emissions")
    a.set_xlabel("x0 [mm]"); a.set_ylabel("y0 [mm]"); a.set_aspect("equal")
    a.set_title("Source fill (transverse)")

    # 7. source fill, axial
    a = ax(7)
    limz = 1.05 * max(np.abs(d["z0_mm"]).max(), 1.0)
    h = a.hist2d(d["x0_mm"], d["z0_mm"], bins=60, range=[[-lim, lim], [-limz, limz]], cmap="magma")
    fig.colorbar(h[3], ax=a, label="emissions")
    a.set_xlabel("x0 [mm]"); a.set_ylabel("z0 [mm]")
    a.set_title("Source fill (axial)")

    # 8. coincidence Δt = t1 − t2 (RAW, not TOF-corrected) — the quantity the τ window cuts on:
    #    trues/scatter peak at 0 (the TOF spread), randoms are flat across the window. The TOF-
    #    corrected residual (dt_ns = jitter only) is the timing resolution, annotated as one number.
    a = ax(8)
    draw = d["t1_ns"] - d["t2_ns"]
    xr = tau if tau > 0 else float(np.nanpercentile(np.abs(draw), 99.5) or 1.0)
    for mask, c, lab in ((is_t, C_TRUE, "true"), (is_s, C_SCAT, "scatter"), (is_r, C_RAND, "random")):
        if mask.any():
            a.hist(draw[mask], bins=70, range=(-xr, xr), histtype="step", color=c, label=lab)
    a.axvline(0, color="gray", ls="--", lw=1)
    for s in (-tau, tau):
        if tau > 0:
            a.axvline(s, color="gray", ls=":", lw=1)
    a.set_yscale("log"); a.set_xlabel("Δt = t1 − t2 [ns]"); a.set_ylabel("LORs"); a.legend()
    a.set_title("Coincidence Δt (raw, τ window)")
    dt = d["dt_ns"]; fin = np.isfinite(dt)
    med = np.median(np.abs(dt[is_t & fin])) if (is_t & fin).any() else float("nan")
    a.text(0.03, 0.97, f"TOF-corr.\nmedian|dt|\n{med:.3f} ns", transform=a.transAxes,
           ha="left", va="top", fontsize=8, family="monospace",
           bbox=dict(boxstyle="round", fc="white", ec="0.7"))

    # 9. composition + efficiency: true / scatter(single+multiple) / random
    a = ax(9)
    n_single = int(((nsc == 1) & is_s).sum())
    n_multi = int(((nsc >= 2) & is_s).sum())
    a.bar("true", n_t, color=C_TRUE)
    a.bar("scatter", n_single, color=C_SINGLE, label="single")
    a.bar("scatter", n_multi, bottom=n_single, color=C_MULTI, label="multiple")
    a.bar("random", n_r, color=C_RAND)
    a.set_ylabel("LORs"); a.legend(loc="upper right", fontsize=8)
    a.set_title("Composition")
    acc = n / nev if nev else float("nan")
    txt = (f"N annihilations = {nev:,}\n"
           f"n LORs = {n:,}\n"
           f"acceptance = {acc:.2%}\n"
           f"true    {n_t/max(n,1):.1%}\n"
           f"scatter {n_s/max(n,1):.1%}  (S {n_single/max(n,1):.1%} / M {n_multi/max(n,1):.1%})\n"
           f"random  {n_r/max(n,1):.1%}\n"
           f"R/(T+S) = {n_r/max(n_t+n_s,1):.2%}\n"
           f"σ_xyz {sig:g} mm")
    a.text(0.97, 0.97, txt, transform=a.transAxes, ha="right", va="top", fontsize=8,
           family="monospace", bbox=dict(boxstyle="round", fc="white", ec="0.7"))

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=120)
    print(f"wrote {args.out}")
    print(f"  {n:,} LORs;  acceptance {acc:.2%};  "
          f"true {n_t/max(n,1):.1%} / scatter {n_s/max(n,1):.1%} / random {n_r/max(n,1):.1%}")


if __name__ == "__main__":
    main()
