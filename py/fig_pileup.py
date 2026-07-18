#!/usr/bin/env python3
"""Figure for latex/dead_time.tex — pileup probability vs time after irradiation end.

Reads studies/pileup/{activity_curve.csv,pileup_table.csv} (written by
scripts/studies/pileup_study.jl) and writes latex/figures/pileup_vs_time.png:
left, the scenario activity A(t) at 1 Gy (total + per isotope); right, the per-LOR
pileup probability on the hottest block, P_LOR(t) = 1 - exp(-2 r_hot(t) T_int) with
r_hot(t) = A(t) * Y * s_hot, one curve per scanner x crystal readout, with the three
300 s acquisition windows shaded.

    python3 py/fig_pileup.py
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STUDY  = os.path.join(REPO, "studies", "pileup")
FIGS   = os.path.join(REPO, "latex", "figures")
T_DEL  = [120.0, 180.0, 300.0]
T_AC   = 300.0

# color by crystal readout (= integration time), linestyle by scanner family
CRYST_COLOR = {"CsI": "#1f77b4", "BGO_195K": "#ff7f0e", "BGO_77K": "#d62728"}
SCAN_STYLE  = {"ring 1 m": "-", "r35 50 cm": "--", "r40 50 cm": "--",
               "r35 35 cm": ":", "r40 35 cm": ":", "CHS": "-."}


def read_rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def main():
    os.makedirs(FIGS, exist_ok=True)
    act  = read_rows(os.path.join(STUDY, "activity_curve.csv"))
    tab  = read_rows(os.path.join(STUDY, "pileup_table.csv"))
    t    = np.array([float(r["t_s"]) for r in act])
    Atot = np.array([float(r["A_total_kBq"]) for r in act])   # kBq at 1 Gy

    fig, (axA, axP) = plt.subplots(1, 2, figsize=(12, 4.6))

    isos = [k[2:-4] for k in act[0] if k.startswith("A_") and k != "A_total_kBq"]
    axA.plot(t, Atot, "k-", lw=2, label="total")
    for iso in isos:
        axA.plot(t, [float(r[f"A_{iso}_kBq"]) for r in act], lw=1.2, label=iso)
    axA.set_yscale("log")
    axA.set_ylim(1e-1, 1e3)
    axA.set_xlabel("t after irradiation end [s]")
    axA.set_ylabel("activity [kBq]")
    axA.set_title("Scenario activity, 1 Gy (fast budget)")
    axA.legend(fontsize=8, ncol=2)

    for r in tab:
        Y, s_hot = float(r["singles_per_decay"]), float(r["hot_frac"])
        T_us     = float(r["T_int_us"])
        p = 100 * (1 - np.exp(-2 * (Atot * 1e3) * Y * s_hot * T_us * 1e-6))
        axP.plot(t, p, SCAN_STYLE[r["label"]], color=CRYST_COLOR[r["crystal"]], lw=1.6,
                 label=f'{r["label"]} {r["crystal"].replace("_", " ")} ({T_us:.0f} µs)')
    for td in T_DEL:
        axP.axvspan(td, td + T_AC, color="gray", alpha=0.06)
    axP.axhline(1.0, color="gray", lw=0.8, ls="--")
    axP.set_yscale("log")
    axP.set_xlabel("t after irradiation end [s]")
    axP.set_ylabel(r"$P_{\rm LOR}$ pileup, hottest block [%]")
    axP.set_title("Per-LOR pileup vs time (acquisition windows shaded)")
    axP.legend(fontsize=7, ncol=2)

    for ax in (axA, axP):
        ax.grid(alpha=0.3, which="both")
        ax.set_xlim(0, 600)
    fig.tight_layout()
    out = os.path.join(FIGS, "pileup_vs_time.png")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
