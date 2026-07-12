#!/usr/bin/env python3
"""Figure + table for latex/scanner_prods.tex — the CsI three-ring v2 survey.

Reads the three CsI v2 masters (crysp_ring_1m / crysp_r35_50cm / crysp_r35_35cm, CsI, tumour-centred,
t_ac = 300 s) across the three acquisition scenarios (del120/180/300) and summarises the geometry
dependence: acceptance (fraction of annihilations that yield a selected LOR) and purity (true
fraction of the detected LORs) vs the axial FOV, per scenario. Uses shard 0 of each leaf (acceptance
and purity are per-shard-stable; the master pools ten such shards). Prints the LaTeX table body and
writes latex/figures/csi_survey.png.

    python3 py/fig_csi_survey.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py

REPO  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRODS = os.environ.get("PRODS", os.path.expanduser("~/Projects/PtCryspProds"))
SCENDIR = os.path.join(PRODS, "uniform_headep_sobp_1e8")
FIGS  = os.path.join(REPO, "latex", "figures")
RINGS = [("crysp_ring_1m_csi_2x0", "R38.7/102"), ("crysp_r35_50cm_csi_2x0", "R35/51"),
         ("crysp_r35_35cm_csi_2x0", "R35/36")]
SCEN  = [("del120s_ac300s_1Gy", 120), ("del180s_ac300s_1Gy", 180), ("del300s_ac300s_1Gy", 300)]
SCEN_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c"]


def _read(scanner, leaf):
    p = os.path.join(SCENDIR, scanner, "csi", leaf, "lors_shard000.h5")
    with h5py.File(p, "r") as f:
        nev = int(f.attrs["nevents"]); nrows = int(f.attrs["nrows"]); afov = float(f.attrs["afov_mm"])
        tr = f["truth"][:]
        nt = int(np.count_nonzero(tr == 0)); ns = int(np.count_nonzero(tr == 1)); nr = int(np.count_nonzero(tr == 2))
    return dict(nev=nev, nrows=nrows, afov=afov, acc=100 * nrows / nev, pur=100 * nt / nrows,
                nt=nt, ns=ns, nr=nr)


def main():
    os.makedirs(FIGS, exist_ok=True)
    data = {(scanner, td): _read(scanner, leaf) for scanner, _ in RINGS for leaf, td in SCEN}

    afovs = [_read(sc, SCEN[0][0])["afov"] / 10.0 for sc, _ in RINGS]   # cm

    fig, (axA, axP) = plt.subplots(1, 2, figsize=(12, 4.6))
    for (leaf, td), col in zip(SCEN, SCEN_COLORS):
        accs = [data[(sc, td)]["acc"] for sc, _ in RINGS]
        purs = [data[(sc, td)]["pur"] for sc, _ in RINGS]  # noqa
        axA.plot(afovs, accs, "o-", color=col, lw=1.8, label=r"$t_{\rm del}=%d$ s" % td)
        axP.plot(afovs, purs, "o-", color=col, lw=1.8, label=r"$t_{\rm del}=%d$ s" % td)
    for ax in (axA, axP):
        ax.set_xlabel("axial FOV [cm]"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
        for x, (_, lab) in zip(afovs, RINGS):
            ax.annotate(lab, (x, ax.get_ylim()[0]), fontsize=7, ha="center", va="bottom", alpha=0.6)
    axA.set_ylabel("acceptance [%]"); axA.set_title("Acceptance vs axial FOV (CsI, per scenario)")
    axP.set_ylabel("true fraction of detected LORs [%]"); axP.set_title("Purity vs axial FOV")
    axP.set_ylim(80, 92)
    fig.suptitle("CsI three-ring v2 survey — tumour-centred, $t_{\\rm ac}=300$ s (shard 0)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = os.path.join(FIGS, "csi_survey.png")
    fig.savefig(out, dpi=120); plt.close(fig); print("wrote", out)

    # LaTeX table body (per shard; ×10 for the master).
    print("\n% --- table body (per shard 0; master = 10 shards) ---")
    for scanner, lab in RINGS:
        for leaf, td in SCEN:
            d = data[(scanner, td)]
            print(r"%s & %d & %.0f & %.2f\,M & %.2f & %.1f \\" %
                  (lab, td, d["afov"] / 10, d["nrows"] / 1e6, d["acc"], d["pur"]))


if __name__ == "__main__":
    main()
