#!/usr/bin/env python3
"""Figure + BGO table for latex/scanner_prods.tex — CsI vs BGO_195K, v2 masters.

Reads the three CsI masters (crysp_ring_1m/r35_50cm/r35_35cm, r 35–38.7 cm) and the three BGO
masters (crysp_ring_1m/r40_50cm/r40_35cm, r 40–43.7 cm = CsI + 5 cm cryostat), tumour-centred,
t_ac = 300 s. The CsI and BGO rings share the same axial FOVs (36/51/102 cm), so acceptance and
purity plot against AFOV with one CsI series and one BGO series per scenario. Uses shard 0 of each
leaf (acceptance/purity are per-shard-stable). Prints the LaTeX table body for the BGO survey and
writes latex/figures/crystal_compare.png.

    python3 py/fig_crystal_compare.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import h5py

REPO  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRODS = os.environ.get("PRODS", os.path.expanduser("~/Projects/PtCryspProds"))
SCENDIR = os.path.join(PRODS, "uniform_headep_sobp_1e8")
FIGS  = os.path.join(REPO, "latex", "figures")
SCEN  = [("del120s_ac300s_1Gy", 120), ("del180s_ac300s_1Gy", 180), ("del300s_ac300s_1Gy", 300)]
SCEN_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c"]
# crystal -> (products subdir, [(scanner, R-label), ...] ordered largest→smallest AFOV)
CRYST = {
    "CsI":      ("csi",      [("crysp_ring_1m_csi_2x0", "R38.7"), ("crysp_r35_50cm_csi_2x0", "R35"),
                              ("crysp_r35_35cm_csi_2x0", "R35")]),
    "BGO_195K": ("bgo_195k", [("crysp_ring_1m_bgo_2x0", "R43.7"), ("crysp_r40_50cm_bgo_2x0", "R40"),
                              ("crysp_r40_35cm_bgo_2x0", "R40")]),
}
STYLE = {"CsI": ("-", "o"), "BGO_195K": ("--", "s")}


def _read(scanner, sub, leaf):
    p = os.path.join(SCENDIR, scanner, sub, leaf, "lors_shard000.h5")
    with h5py.File(p, "r") as f:
        nev = int(f.attrs["nevents"]); nrows = int(f.attrs["nrows"]); afov = float(f.attrs["afov_mm"]) / 10.0
        nt = int(np.count_nonzero(f["truth"][:] == 0))
    return dict(acc=100 * nrows / nev, pur=100 * nt / nrows, nrows=nrows, afov=afov)


def main():
    os.makedirs(FIGS, exist_ok=True)
    D = {(cr, ring, td): _read(scanner, sub, leaf)
         for cr, (sub, rings) in CRYST.items()
         for ring, (scanner, _) in enumerate(rings)
         for leaf, td in SCEN}

    fig, (axA, axP) = plt.subplots(1, 2, figsize=(12, 4.8))
    for cr, (ls, mk) in STYLE.items():
        afovs = [D[(cr, r, 120)]["afov"] for r in range(3)]
        for (leaf, td), col in zip(SCEN, SCEN_COLORS):
            axA.plot(afovs, [D[(cr, r, td)]["acc"] for r in range(3)], ls, marker=mk, color=col, lw=1.7, ms=6)
            axP.plot(afovs, [D[(cr, r, td)]["pur"] for r in range(3)], ls, marker=mk, color=col, lw=1.7, ms=6)
    axA.set_ylabel("acceptance [%]"); axA.set_title("Acceptance vs axial FOV — CsI vs BGO")
    axP.set_ylabel("true fraction of detected LORs [%]"); axP.set_title("Purity vs axial FOV")
    axP.set_ylim(68, 92)
    # legends: colour = scenario, style = crystal
    scen_h = [mlines.Line2D([], [], color=c, lw=2, label=r"$t_{\rm del}=%d$ s" % td) for (_, td), c in zip(SCEN, SCEN_COLORS)]
    cryst_h = [mlines.Line2D([], [], color="0.35", ls=ls, marker=mk, label=cr) for cr, (ls, mk) in STYLE.items()]
    axA.legend(handles=scen_h, fontsize=8, loc="upper left")
    axP.legend(handles=cryst_h, fontsize=8, loc="upper right")
    for ax in (axA, axP):
        ax.set_xlabel("axial FOV [cm]"); ax.grid(alpha=0.3)
    fig.suptitle("CsI vs BGO$_{195\\rm K}$ v2 masters — tumour-centred, $t_{\\rm ac}=300$ s "
                 "(BGO at $R_{\\rm CsI}+5$ cm cryostat; shard 0)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(FIGS, "crystal_compare.png")
    fig.savefig(out, dpi=120); plt.close(fig); print("wrote", out)

    print("\n% --- BGO survey table body (per shard 0; master = 10 shards) ---")
    labels = ["R43.7/102", "R40/51", "R40/36"]
    for ring, lab in enumerate(labels):
        for leaf, td in SCEN:
            d = D[("BGO_195K", ring, td)]
            print(r"%s & %d & %.0f & %.2f\,M & %.2f & %.1f \\" % (lab, td, d["afov"], d["nrows"] / 1e6, d["acc"], d["pur"]))


if __name__ == "__main__":
    main()
