#!/usr/bin/env python3
"""Source-level truth figures for latex/scanner_prods.tex (v2 convention).

Both figures are SOURCE-level (truth, no scanner), in the v2 convention: **tumour-centred** — z = 0
is the SOBP dose-target centre (a fixed anatomical reference, from the distal dose R80 minus half the
target thickness), not the drifting activity edge. The acquisition is a family of **fixed windows**
[t_del, t_del + t_ac] on the irradiation-end clock (t_ac = 300 s), one per patient-arrival delay
t_del. The activity a scanner would record in scenario (t_del) is the truth profile with each isotope
j scaled by its expected decays in that window,

    s_j(t1, t2) = (e^{-lam_j t1} - e^{-lam_j t2}) / (e^{-lam_j t_del^bud} - e^{-lam_j (t_del^bud + t_meas^bud)}),

i.e. the count in [t1, t2] relative to the fast budget window the truth profile is normalised to
([60, 1260] s). Uses only the truth activity profile + dose + the isotope half-lives — no scanner.

Produces two tracked figures in latex/figures/:
  activity_truth.png       the truth activity per isotope (full window), tumour centred at 0.
  activity_vs_tstart.png   three panels, one per fixed-window scenario t_del = 120/180/300 s.

Source: <PRODS>/<scenario>/truth/{activity_profile_<budget>.csv, depth_dose.csv, run_meta.csv}.

    python3 py/fig_activity_profiles.py
"""
import csv
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRODS  = os.environ.get("PRODS", os.path.expanduser("~/Projects/PtCryspProds"))
SCEN   = "uniform_headep_sobp_1e8"
BUDGET = "fast"
FIGS   = os.path.join(REPO, "latex", "figures")
TRUTH  = os.path.join(PRODS, SCEN, "truth")

# Isotope half-lives [s] — frozen scenario constants (isotopes.csv).
HALF_LIVES = {"O15": 122.24, "C11": 1223.4, "N13": 597.9, "C10": 19.29, "O14": 70.62}
# The fast budget window [t_del_bud, t_del_bud + t_meas_bud] the truth profile counts are normalised to.
T_DEL_BUD, T_MEAS_BUD = 60.0, 1200.0
# v2 acquisition scenarios: fixed length t_ac, patient-arrival delays t_del (irradiation-end clock).
T_AC   = 300.0
TDELS  = [120, 180, 300]
ISO_COLORS = {"O15": "#226688", "C11": "#ee88aa", "N13": "#44aa44", "C10": "#cc3333", "O14": "#aa7700"}


def _window_scale(iso, t1, t2):
    """Expected decays of `iso` in [t1, t2] relative to the fast budget window (the profile's norm)."""
    lam = np.log(2) / HALF_LIVES[iso]
    return (np.exp(-lam * t1) - np.exp(-lam * t2)) / (np.exp(-lam * T_DEL_BUD) - np.exp(-lam * (T_DEL_BUD + T_MEAS_BUD)))


def _tumour_center_mm():
    """Native-frame z [mm] of the SOBP target centre = distal dose R80 − half the target thickness."""
    rows = list(csv.DictReader(open(os.path.join(TRUTH, "depth_dose.csv"))))
    z = np.array([float(r["z_mm"]) for r in rows]); d = np.array([float(r["dose_core_Gy"]) for r in rows])
    im = int(np.argmax(d)); half = 0.8 * d[im]; zd80 = z[im]
    for i in range(im, len(z) - 1):
        if d[i] >= half > d[i + 1]:
            zd80 = z[i] + (z[i + 1] - z[i]) * (d[i] - half) / (d[i] - d[i + 1]); break
    m = list(csv.DictReader(open(os.path.join(TRUTH, "run_meta.csv"))))[0]
    prox, dist = float(m["target_prox_depth_mm"]), float(m["target_dist_depth_mm"])
    return zd80 - 0.5 * (dist - prox), 0.5 * (dist - prox)     # (centre, half-thickness)


def _load_profile():
    """Return (z tumour-centred, {iso: activity(z)}, isos, tumour half-thickness [mm])."""
    rows = list(csv.DictReader(open(os.path.join(TRUTH, "activity_profile_%s.csv" % BUDGET))))
    isos = [c for c in rows[0] if c not in ("z_mm", "total")]
    z    = np.array([float(r["z_mm"]) for r in rows])
    A    = {iso: np.array([float(r[iso]) for r in rows]) for iso in isos}
    center, halfthick = _tumour_center_mm()
    return z - center, A, isos, halfthick                     # tumour → z = 0


def fig_truth():
    """activity_truth.png — the truth activity per isotope (full window), tumour centred at 0."""
    z, A, isos, ht = _load_profile()
    fig, ax = plt.subplots(figsize=(8, 4.6))
    for iso in isos:
        ax.plot(z, A[iso], lw=1.3, color=ISO_COLORS.get(iso), label=iso)
    ax.plot(z, sum(A.values()), lw=2.2, color="k", label="total")
    ax.axvspan(-ht, ht, alpha=0.12, color="purple", label="tumour/target extent")
    ax.axvline(0, ls="--", color="k", lw=1, label="tumour centre = ring centre")
    ax.set_xlabel("z along beam [mm]  (tumour/target centre at ring centre)")
    ax.set_ylabel("expected decays / bin")
    ax.set_title("Source activity depth profile (truth), per isotope — 1 Gy")
    ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.3); ax.set_xlim(-90, 40)
    fig.tight_layout(); out = os.path.join(FIGS, "activity_truth.png")
    fig.savefig(out, dpi=110); plt.close(fig); print("wrote", out)


def fig_tstart():
    """activity_vs_tstart.png — three panels, per fixed-window scenario t_del = 120/180/300 s."""
    z, A, isos, ht = _load_profile()
    ymax = max(sum(A.values()))
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharex=True, sharey=True)
    for ax, td in zip(axes, TDELS):
        t1, t2 = float(td), float(td) + T_AC
        scaled = {iso: A[iso] * _window_scale(iso, t1, t2) for iso in isos}
        for iso in isos:
            ax.plot(z, scaled[iso], lw=1.3, color=ISO_COLORS.get(iso), label=iso)
        ax.plot(z, sum(scaled.values()), lw=2.2, color="k", label="total")
        ax.axvline(0, ls=":", color="k", lw=0.9)
        ax.set_title(r"$t_{\rm del}=%d$ s, window $[%d,%d]$ s" % (td, td, td + int(T_AC)), fontsize=10)
        ax.grid(alpha=0.3); ax.set_xlim(-90, 40); ax.set_ylim(0, 1.05 * ymax)
        ax.set_xlabel("z along beam [mm]  (tumour at 0)")
    axes[0].set_ylabel("expected decays / bin"); axes[0].legend(fontsize=8, ncol=2)
    fig.suptitle("Source activity per fixed acquisition scenario (truth) — 1 Gy, $t_{\\rm ac}=300$ s", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96]); out = os.path.join(FIGS, "activity_vs_tstart.png")
    fig.savefig(out, dpi=110); plt.close(fig); print("wrote", out)


def main():
    os.makedirs(FIGS, exist_ok=True)
    fig_truth()
    fig_tstart()


if __name__ == "__main__":
    main()
