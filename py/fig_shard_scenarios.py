#!/usr/bin/env python3
"""Figure for latex/scanner_prods.tex — the first v2 shard round (DETECTED LORs, per scenario).

Reads the three acquisition-scenario masters of the reference detector
(crysp_ring_1m_csi_2x0 / csi / del{120,180,300}s_ac300s_1Gy) and shows what the scanner actually
records for each patient-arrival delay, in the v2 convention: tumour-centred (z = 0 is the SOBP
dose-target centre) and the emitting isotope carried per LOR. Two panels:

  (left)  detected activity depth profile a(z0) per scenario (trues only), z0 = the true
          annihilation depth, already tumour-centred in the shard. Later start → fewer counts and a
          reshaped distal edge.
  (right) isotope composition of the accepted trues per scenario (the `isotope` column): O15 depletes
          and C11 grows as the start is delayed — the isotope-mix evolution, now measurable per LOR.

Source: <PRODS>/uniform_headep_sobp_1e8/crysp_ring_1m_csi_2x0/csi/<leaf>/lors_shard*.h5 (all shards
in each leaf are pooled). Writes latex/figures/shard_scenarios.png.

    python3 py/fig_shard_scenarios.py
"""
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py

REPO  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRODS = os.environ.get("PRODS", os.path.expanduser("~/Projects/PtCryspProds"))
LEAF  = os.path.join(PRODS, "uniform_headep_sobp_1e8", "crysp_ring_1m_csi_2x0", "csi")
FIGS  = os.path.join(REPO, "latex", "figures")
SCEN  = [("del120s_ac300s_1Gy", 120), ("del180s_ac300s_1Gy", 180), ("del300s_ac300s_1Gy", 300)]
ISO   = ["O15", "C11", "N13", "C10", "O14"]
ISO_COLORS = {"O15": "#226688", "C11": "#ee88aa", "N13": "#44aa44", "C10": "#cc3333", "O14": "#aa7700"}
SCEN_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c"]
ZEDGES = np.arange(-90.0, 40.01, 2.0)                      # tumour-centred z bins [mm]


def _pool(leaf):
    """Pool all shards in a leaf → (z0 tumour-centred [mm], isotope) for trues only, + attrs."""
    files = sorted(glob.glob(os.path.join(LEAF, leaf, "lors_shard*.h5")))
    if not files:
        raise SystemExit("no shards in %s — run the v2 master first" % os.path.join(LEAF, leaf))
    zs, isos = [], []
    attrs = {}
    for p in files:
        with h5py.File(p, "r") as f:
            if not attrs:
                attrs = {k: f.attrs[k] for k in ("source_z_offset_mm", "washout_g", "t_del_s", "t_ac_s", "center_on")}
            m = f["truth"][:] == 0                          # trues only
            zs.append(f["z0_mm"][:][m].astype(np.float64) * float(f.attrs["xyz_scale_mm"]))
            isos.append(f["isotope"][:][m])
    return np.concatenate(zs), np.concatenate(isos), attrs


def main():
    os.makedirs(FIGS, exist_ok=True)
    data = [( _pool(leaf), td) for leaf, td in SCEN]

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.5, 4.8))

    # ---- left: detected depth profile per scenario (tumour-centred) ----
    zc = 0.5 * (ZEDGES[:-1] + ZEDGES[1:])
    for ((z, iso, at), td), col in zip(data, SCEN_COLORS):
        h, _ = np.histogram(z, bins=ZEDGES)
        axL.step(zc, h, where="mid", lw=1.6, color=col,
                 label=r"$t_{\rm del}=%d$ s  (%.2fM trues)" % (td, len(z) / 1e6))
    axL.axvline(0, ls="--", color="k", lw=1.1, label="tumour centre ($z=0$)")
    axL.set_xlabel("z along beam [mm]  (tumour/target centre at 0)")
    axL.set_ylabel("detected trues / 2 mm bin")
    axL.set_title("Detected activity depth profile, per acquisition scenario")
    axL.legend(fontsize=8); axL.grid(alpha=0.3); axL.set_xlim(ZEDGES[0], ZEDGES[-1])

    # ---- right: isotope composition of the detected trues per scenario ----
    x = np.arange(len(SCEN)); w = 0.8
    fracs = np.zeros((len(SCEN), len(ISO)))
    for i, ((z, iso, at), td) in enumerate(data):
        cnt = np.array([np.count_nonzero(iso == j) for j in range(len(ISO))], float)
        fracs[i] = cnt / max(cnt.sum(), 1)
    bottom = np.zeros(len(SCEN))
    for j, name in enumerate(ISO):
        axR.bar(x, fracs[:, j], w, bottom=bottom, color=ISO_COLORS[name], label=name)
        bottom += fracs[:, j]
    axR.set_xticks(x); axR.set_xticklabels([r"$t_{\rm del}=%d$ s" % td for _, td in SCEN])
    axR.set_ylabel("fraction of detected trues")
    axR.set_title("Isotope composition (from the per-LOR isotope truth)")
    axR.set_ylim(0, 1); axR.legend(fontsize=8, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.10))
    axR.grid(alpha=0.3, axis="y")

    off = float(data[0][0][2]["source_z_offset_mm"])
    fig.suptitle("First v2 master round — crysp_ring_1m CsI, tumour-centred (offset +%.1f mm), $t_{\\rm ac}=300$ s"
                 % off, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(FIGS, "shard_scenarios.png")
    fig.savefig(out, dpi=120); plt.close(fig); print("wrote", out)


if __name__ == "__main__":
    main()
