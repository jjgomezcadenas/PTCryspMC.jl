#!/usr/bin/env python3
"""Plot the activity-model event-time histogram (scripts/studies/o15_lifetime.jl output)
against the analytic truncated-exponential decay curve.

    python3 py/plot_o15_lifetime.py                       # output/o15_lifetime.csv
    python3 py/plot_o15_lifetime.py --csv path.csv --out plot.png
"""
import argparse
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path):
    """Return (params dict, t_lo, t_hi, count) — params from the leading '# k=v' header, the
    rest a single manual pass (robust to the comment + named-header combo)."""
    params, rows = {}, []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                for tok in s[1:].split():
                    if "=" in tok:
                        k, v = tok.split("=", 1)
                        try:
                            params[k] = float(v)
                        except ValueError:
                            params[k] = v
            elif s.startswith("t_lo"):
                continue  # the column-name header row
            else:
                a, b, c = s.split(",")
                rows.append((float(a), float(b), float(c)))
    arr = np.array(rows)
    return params, arr[:, 0], arr[:, 1], arr[:, 2]


def main():
    ap = argparse.ArgumentParser(description="Plot the activity-model event-time distribution.")
    ap.add_argument("--csv", default="output/o15_lifetime.csv")
    ap.add_argument("--out", default="")
    a = ap.parse_args()

    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    csv = a.csv if os.path.isabs(a.csv) else os.path.join(repo, a.csv)
    out = a.out or os.path.splitext(csv)[0] + ".png"
    out = out if os.path.isabs(out) else os.path.join(repo, out)

    params, t_lo, t_hi, count = read_csv(csv)
    lam = params["lambda_per_s"]; t0 = params["t0_s"]; t1 = params["t1_s"]
    half = params.get("half_life_s", math.log(2) / lam)
    N = count.sum(); w = t_hi[0] - t_lo[0]
    centers = 0.5 * (t_lo + t_hi)

    # Analytic truncated-exponential density on [t0,t1] (rate lam), scaled to counts-per-bin.
    norm = 1.0 - math.exp(-lam * (t1 - t0))
    tt = np.linspace(t0, t1, 400)
    pdf = lam * np.exp(-lam * (tt - t0)) / norm           # 1/s
    expected_per_bin = pdf * w * N

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for ax, logy in zip(axes, (False, True)):
        ax.bar(centers, count, width=w, align="center", color="#4C72B0", alpha=0.65,
               label="sampled events")
        ax.plot(tt, expected_per_bin, color="#C44E52", lw=2, label="analytic $\\propto e^{-\\lambda t}$")
        ax.axvline(t0 + half, color="gray", ls="--", lw=1, label=f"T½ = {half:.0f} s")
        ax.set_xlabel("annihilation time [s]")
        ax.set_ylabel("events / bin")
        if logy:
            ax.set_yscale("log"); ax.set_title("log scale (straight line = exponential)")
        else:
            ax.set_title("linear scale (front-loaded)")
        ax.legend(frameon=False)

    fig.suptitle(f"Activity-model event times — ¹⁵O-style T½={half:.1f} s, window [{t0:.0f}, {t1:.0f}] s, "
                 f"N={int(N):,}", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
