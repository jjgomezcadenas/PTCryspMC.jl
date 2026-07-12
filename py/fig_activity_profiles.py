#!/usr/bin/env python3
"""Figures for latex/scanner_prods.tex — the source activity and its acquisition-time dependence.

Both figures are SOURCE-level (truth, no scanner): the beta+ activity depth profile per isotope,
with the distal edge (R50) centred at z=0 (patient positioning). At acquisition start time
t_start the activity a scanner records is the truth profile with each isotope j scaled by its
surviving fraction over the window [0, T_meas],

    f_j(t_start) = (e^{-lam_j t_start} - e^{-lam_j T_meas}) / (1 - e^{-lam_j T_meas}),

since the decay time is drawn from a truncated exponential on [0, T_meas] independent of position
(so a start-time cut scales each isotope's amplitude, not its shape). This uses only the truth
activity profile and the isotope half-lives — no per-scanner simulation.

Produces two tracked figures in latex/figures/:
  activity_truth.png       the truth activity per isotope (full window), edge centred at 0.
  activity_vs_tstart.png   four panels, one per t_start = 0,120,180,300 s, each the per-isotope
                           activity scaled by f_j(t_start) — the isotope-composition evolution.

Source: <PRODS>/<scenario>/truth/activity_profile_<budget>.csv (the canonical truth product from
src/scenario.jl:write_activity_profile). Half-lives / window from the scenario (below).

    python3 py/fig_activity_profiles.py
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRODS = os.environ.get("PRODS", os.path.expanduser("~/Projects/PtCryspProds"))
SCEN  = "uniform_headep_sobp_1e8"
BUDGET = "fast"
FIGS  = os.path.join(REPO, "latex", "figures")

# Isotope half-lives [s] and the acquisition window [s] — frozen scenario constants
# (isotopes.csv / sampling_budget_<budget>_meta.csv of the scenario).
HALF_LIVES = {"O15": 122.24, "C11": 1223.4, "N13": 597.9, "C10": 19.29, "O14": 70.62}
T_MEAS = 1200.0
TSTARTS = [0, 120, 180, 300]
ISO_COLORS = {"O15": "#226688", "C11": "#ee88aa", "N13": "#44aa44", "C10": "#cc3333", "O14": "#aa7700"}


def _distal_edge(z, y):
    """Downstream (+z) 50%-of-peak crossing of profile y(z) — the activity R50 [mm]."""
    pk = int(np.argmax(y)); half = y[pk] / 2.0
    for i in range(pk, len(y) - 1):
        if y[i] >= half > y[i + 1]:
            return z[i] + (z[i + 1] - z[i]) * (y[i] - half) / (y[i] - y[i + 1])
    return z[pk]


def _surviving(iso, t_start):
    """Fraction of isotope `iso` decaying in [t_start, T_MEAS] (truncated-exponential window)."""
    lam = np.log(2) / HALF_LIVES[iso]
    return (np.exp(-lam * t_start) - np.exp(-lam * T_MEAS)) / (1.0 - np.exp(-lam * T_MEAS))


# Tumour/target z-extent in the native scenario frame [mm] (headep tumour ellipsoid, c=20 at -25).
TUMOUR_MM = (-45.0, -5.0)


def _load_profile():
    """Return (z centred on the distal edge, {iso: activity(z)}, isos, shift[mm]) from the bundle."""
    path = os.path.join(PRODS, SCEN, "truth", "activity_profile_%s.csv" % BUDGET)
    rows = list(csv.DictReader(open(path)))
    isos = [c for c in rows[0] if c not in ("z_mm", "total")]
    z    = np.array([float(r["z_mm"]) for r in rows])
    A    = {iso: np.array([float(r[iso]) for r in rows]) for iso in isos}
    tot  = np.array([float(r["total"]) for r in rows])
    shift = -_distal_edge(z, tot)                      # bring the distal edge to z=0
    return z + shift, A, isos, shift


def fig_truth():
    """activity_truth.png — the truth activity per isotope (full window), edge centred at 0."""
    z, A, isos, shift = _load_profile()
    fig, ax = plt.subplots(figsize=(8, 4.6))
    for iso in isos:
        ax.plot(z, A[iso], lw=1.3, color=ISO_COLORS.get(iso), label=iso)
    ax.plot(z, sum(A.values()), lw=2.2, color="k", label="total")
    ax.axvspan(TUMOUR_MM[0] + shift, TUMOUR_MM[1] + shift, alpha=0.12, color="purple", label="tumour region")
    ax.axvline(0, ls="--", color="k", lw=1, label="distal edge (R50) = ring centre")
    ax.set_xlabel("z along beam [mm]  (distal edge at ring centre)")
    ax.set_ylabel("expected decays / bin")
    ax.set_title("Source activity depth profile (truth), per isotope — 1 Gy")
    ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3); ax.set_xlim(-100, 25)
    fig.tight_layout(); out = os.path.join(FIGS, "activity_truth.png")
    fig.savefig(out, dpi=110); plt.close(fig); print("wrote", out)


def fig_tstart():
    """activity_vs_tstart.png — four panels (t_start=0,120,180,300 s), per-isotope source activity."""
    z, A, isos, _ = _load_profile()
    ymax = max(sum(A.values()))
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True, sharey=True)
    axes = axes.flatten()
    for ax, t in zip(axes, TSTARTS):
        scaled = {iso: A[iso] * _surviving(iso, t) for iso in isos}
        for iso in isos:
            ax.plot(z, scaled[iso], lw=1.3, color=ISO_COLORS.get(iso), label=iso)
        ax.plot(z, sum(scaled.values()), lw=2.2, color="k", label="total")
        ax.axvline(0, ls=":", color="k", lw=0.9)
        ax.set_title("acquisition start $t = %d$ s" % t, fontsize=10)
        ax.grid(alpha=0.3); ax.set_xlim(-100, 25); ax.set_ylim(0, 1.05 * ymax)
    axes[0].legend(fontsize=8, ncol=2)
    for i in (0, 2): axes[i].set_ylabel("expected decays / bin")
    for i in (2, 3): axes[i].set_xlabel("z along beam [mm]  (distal edge at 0)")
    fig.suptitle("Source activity depth profile per isotope, vs acquisition start time — 1 Gy", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.98]); out = os.path.join(FIGS, "activity_vs_tstart.png")
    fig.savefig(out, dpi=110); plt.close(fig); print("wrote", out)


def main():
    os.makedirs(FIGS, exist_ok=True)
    fig_truth()
    fig_tstart()


if __name__ == "__main__":
    main()
