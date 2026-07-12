#!/usr/bin/env python3
"""Figures for latex/scanner_prods.tex — the source activity and its acquisition-time dependence.

Produces two tracked figures in latex/figures/:

  activity_truth.png      the truth beta+ activity depth profile, per isotope, with the distal
                          edge (R50) centred at z=0 (patient positioning). Source:
                          <PRODS>/<scenario>/truth/activity_profile_<budget>.csv (the canonical
                          truth product written by src/scenario.jl:write_activity_profile).

  activity_vs_tstart.png  per scanner, the accepted-LOR activity depth profile at acquisition
                          start times t_start = 0,120,180,300 s (the cut t_decay >= t_start on
                          the stored t_decay_s column). Source: the distal-edge-centred CsI
                          shard-0 lors_det.h5 of each geometry, under <PROD>/<tag>/.

Data prerequisites (regenerable; the lors_det files are gitignored and may be pruned). Each arm
is one config with [source].center_on = "distal_edge", chain simulate -> build_true ->
build_randoms -> reco_lors, realization 0:
    runs/uniform_headep_{ring1m,r30_50cm,ring50cm,r35_50cm,chs,r35_35cm,r30_30cm,r35_30cm}_csi_distal_centered.toml

    python3 py/fig_activity_profiles.py
"""
import csv
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py

REPO  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRODS = os.environ.get("PRODS", os.path.expanduser("~/Projects/PtCryspProds"))
PROD  = os.path.join(REPO, "prod")                 # local production tree (gitignored)
SCEN  = "uniform_headep_sobp_1e8"
BUDGET = "fast"
FIGS  = os.path.join(REPO, "latex", "figures")

# scanner tag -> (prod dir under PROD, display title), in a sensible radius/AFOV order
SCANNERS = [
    ("ring1m",   "uniform_headep_ring1m_csi_distal_centered",   "ring 1m (R387x1024)"),
    ("r30_50cm", "uniform_headep_r30_50cm_csi_distal_centered", "R300x50cm"),
    ("ring50cm", "uniform_headep_ring50cm_csi_distal_centered", "R387x50cm"),
    ("r35_50cm", "uniform_headep_r35_50cm_csi_distal_centered", "R350x50cm"),
    ("chs",      "uniform_headep_chs_csi_distal_centered",      "CHS (R200x358)"),
    ("r35_35cm", "uniform_headep_r35_35cm_csi_distal_centered", "R350x35cm"),
    ("r30_30cm", "uniform_headep_r30_30cm_csi_distal_centered", "R300x30cm"),
    ("r35_30cm", "uniform_headep_r35_30cm_csi_distal_centered", "R350x30cm"),
]
TSTARTS = [(0, "0 s", "#222222"), (120, "120 s", "#226688"),
           (180, "180 s", "#22aa88"), (300, "300 s", "#cc3333")]
ISO_COLORS = {"O15": "#226688", "C11": "#ee88aa", "N13": "#44aa44", "C10": "#cc3333", "O14": "#aa7700"}


def _distal_edge(z, y):
    """Downstream (+z) 50%-of-peak crossing of profile y(z) — the activity R50 [mm]."""
    pk = int(np.argmax(y)); half = y[pk] / 2.0
    for i in range(pk, len(y) - 1):
        if y[i] >= half > y[i + 1]:
            return z[i] + (z[i + 1] - z[i]) * (y[i] - half) / (y[i] - y[i + 1])
    return z[pk]


def fig_truth():
    """activity_truth.png — the truth activity profile, per isotope, centred on its distal edge."""
    path = os.path.join(PRODS, SCEN, "truth", "activity_profile_%s.csv" % BUDGET)
    rows = list(csv.DictReader(open(path)))
    isos = [c for c in rows[0] if c not in ("z_mm", "total")]
    z    = np.array([float(r["z_mm"]) for r in rows])
    tot  = np.array([float(r["total"]) for r in rows])
    shift = -_distal_edge(z, tot)                      # centre the distal edge at z=0
    z = z + shift

    fig, ax = plt.subplots(figsize=(8, 4.6))
    for iso in isos:
        ax.plot(z, [float(r[iso]) for r in rows], lw=1.3, color=ISO_COLORS.get(iso), label=iso)
    ax.plot(z, tot, lw=2.2, color="k", label="total")
    ax.axvspan(-45 + shift, -5 + shift, alpha=0.12, color="purple", label="tumour region")
    ax.axvline(0, ls="--", color="k", lw=1, label="distal edge (R50) = ring centre")
    ax.set_xlabel("z along beam [mm]  (distal edge centred at ring centre)")
    ax.set_ylabel("expected decays / bin")
    ax.set_title("Source activity depth profile (truth), per isotope — 1 Gy")
    ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3); ax.set_xlim(-100, 25)
    fig.tight_layout(); out = os.path.join(FIGS, "activity_truth.png")
    fig.savefig(out, dpi=110); plt.close(fig); print("wrote", out)


def fig_tstart():
    """activity_vs_tstart.png — per scanner, accepted-LOR activity at t_start = 0,120,180,300 s."""
    lo, hi, bw = -100.0, 30.0, 2.0
    edges = np.arange(lo, hi + bw, bw); centers = (edges[:-1] + edges[1:]) / 2
    panels = []
    for tag, d, title in SCANNERS:
        p = os.path.join(PROD, d, "lors_det.h5")
        if not os.path.isfile(p):
            print("  skip (missing):", p); continue
        with h5py.File(p, "r") as f:
            sc = float(f.attrs["xyz_scale_mm"])
            z0 = f["z0_mm"][:].astype(np.float64) * sc
            td = f["t_decay_s"][:].astype(np.float64)
        hists = [np.histogram(z0[td >= t], bins=edges)[0] for t, _, _ in TSTARTS]
        panels.append((title, hists))

    n = len(panels); ncol = 2; nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 2.3 * nrow), sharex=True)
    axes = np.atleast_1d(axes).flatten()
    for ax, (title, hists) in zip(axes, panels):
        for h, (_, lab, c) in zip(hists, TSTARTS):
            ax.plot(centers, h / 1e3, lw=1.4, color=c, label=lab)
        ax.axvline(0, ls=":", color="k", lw=0.8)
        ax.set_title(title, fontsize=9); ax.grid(alpha=0.3); ax.set_xlim(-100, 25)
    for ax in axes[n:]:
        ax.axis("off")
    axes[0].legend(fontsize=7, title="acq. start", title_fontsize=7)
    for i in range(0, n, ncol):
        axes[i].set_ylabel("acc. LORs /2mm [k]")
    for i in range(max(0, n - ncol), n):
        axes[i].set_xlabel("z [mm] (distal edge at 0)")
    fig.suptitle("Accepted-LOR activity vs acquisition start time (CsI, distal-edge centred)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.98]); out = os.path.join(FIGS, "activity_vs_tstart.png")
    fig.savefig(out, dpi=110); plt.close(fig); print("wrote", out, "(%d panels)" % n)


def main():
    os.makedirs(FIGS, exist_ok=True)
    fig_truth()
    fig_tstart()


if __name__ == "__main__":
    main()
