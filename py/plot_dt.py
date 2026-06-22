#!/usr/bin/env python3
"""Plot the timing of the true coincidences (scripts/studies/examine_dt.jl output) to choose the
coincidence window τ.

    dT_ns = |t1 - t2|         the RAW timestamp difference — the coincidence-window variable
    dt_ns = |Δt0| - TOF_diff  the TOF-corrected residual = the timing resolution (diagnostic)

Panels: (1) the dT distribution, (2) the cumulative trues-kept vs τ (read τ off this for a target
acceptance — wider τ keeps more trues but admits more randoms ∝ τ), (3) the dt residual.

    python3 py/plot_dt.py --csv studies/dt/sphere_water_csi_dt.csv
    python3 py/plot_dt.py --csv ... --out plot.png --dtmax 10
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(description="Plot dT/dt of the true coincidences for the τ choice.")
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--dtmax", type=float, default=10.0, help="dT axis upper limit [ns] for the zoomed view")
    a = ap.parse_args()

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    csv = a.csv if os.path.isabs(a.csv) else os.path.join(repo, a.csv)
    out = a.out or os.path.splitext(csv)[0] + ".png"
    out = out if os.path.isabs(out) else os.path.join(repo, out)

    data = np.loadtxt(csv, delimiter=",", skiprows=1)
    dT, dt = data[:, 0], data[:, 1]
    n = dT.size
    qs = [50, 90, 99, 99.9]
    pT = np.percentile(dT, qs)
    print(f"{n:,} LORs   dT_ns percentiles: " + "  ".join(f"p{q}={v:.3f}" for q, v in zip(qs, pT)))

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    # (1) dT distribution, zoomed to the window-relevant region, log-y.
    ax = axes[0]
    ax.hist(dT, bins=200, range=(0, a.dtmax), color="#4C72B0")   # range drops the tail beyond dtmax
    for q, v in zip(qs[2:], pT[2:]):
        if v <= a.dtmax:
            ax.axvline(v, color="#C44E52", ls="--", lw=1, label=f"p{q}={v:.2f} ns")
    ax.set_yscale("log")
    ax.set_xlabel("dT = |t1 - t2|  [ns]")
    ax.set_ylabel("LORs / bin")
    ax.set_title(f"coincidence-time difference (≤ {a.dtmax:g} ns)")
    ax.legend(frameon=False)

    # (2) Cumulative trues kept vs τ — the tradeoff curve.
    ax = axes[1]
    s = np.sort(dT)
    frac = np.arange(1, n + 1) / n
    m = s <= a.dtmax
    ax.plot(s[m], 100 * frac[m], color="#55A868", lw=2)
    for q, v in zip(qs[2:], pT[2:]):
        if v <= a.dtmax:
            ax.axvline(v, color="#C44E52", ls="--", lw=1)
            ax.axhline(q, color="#C44E52", ls=":", lw=0.8)
    ax.set_xlabel("τ  [ns]")
    ax.set_ylabel("trues kept [%]")
    ax.set_ylim(0, 100)
    ax.set_title("acceptance vs window τ")

    # (3) The dt residual = timing resolution.
    ax = axes[2]
    lo, hi = np.percentile(dt, [0.5, 99.5])
    ax.hist(dt, bins=200, range=(lo, hi), color="#8172B3")       # range drops the tails
    ax.set_xlabel("dt = |Δt0| - TOF_diff  [ns]")
    ax.set_ylabel("LORs / bin")
    ax.set_title(f"timing residual  (std {dt.std():.3f} ns)")

    fig.suptitle(f"True-coincidence timing — {os.path.basename(csv)}  (N={n:,})", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
