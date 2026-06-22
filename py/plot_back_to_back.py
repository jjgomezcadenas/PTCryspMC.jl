#!/usr/bin/env python3
"""Summary plots for the back-to-back stack written by
scripts/shoot_into_ring.jl.

Reads a stack CSV (one row per interaction:
event_number, gamma, step, x_mm, y_mm, z_mm, e_in_keV, e_dep_keV, process, iz, iphi)
and draws a 3x3 panel. Quantities are in the crystal-local frame: the depth of
interaction (DOI = r - r_inner, radial, 0..wall) and the crystal face (x = arc
offset across the φ-sector, y = axial offset within the wheel).

  1. Edep per gamma            5. both-γ-contained fraction     9. ring hit map (φ vs z)
  2. 3-D impacts on the ring   6. Edep, contained gammas
  3. DOI of 1st interaction    7. impact x-y within a crystal
  4. DOI of 2nd interaction    8. blocks touched per gamma

"Contained" = a gamma deposits its full 511 keV in a single (iz, iphi) crystal
(no overspill, no escape). The impact = the first energy-depositing interaction.

Run from the repo root:
    python3 py/plot_back_to_back.py --csv output/b2b_csi_stack.csv \
                                    --out output/control_plots/b2b_csi_summary.png
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")  # headless: write a file, never open a window
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KEY = ["event_number", "gamma"]
FULL_KEV = 505.0     # Etot above this counts as the full 511 keV (contained)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", default=os.path.join(REPO, "studies", "b2b", "b2b_csi_stack.csv"))
    p.add_argument("--out", default=os.path.join(REPO, "studies", "b2b",
                                                 "b2b_csi_summary.png"))
    # CRYSP1M ring geometry (for the DOI / crystal-local frame)
    p.add_argument("--r-inner-mm", type=float, default=387.0)
    p.add_argument("--wall-mm", type=float, default=37.0)
    p.add_argument("--half-length-mm", type=float, default=512.0)
    p.add_argument("--n-phi", type=int, default=48)
    p.add_argument("--n-z", type=int, default=20)
    args = p.parse_args()

    Ri, wall, H = args.r_inner_mm, args.wall_mm, args.half_length_mm
    dphi, dz = 2 * np.pi / args.n_phi, 2 * H / args.n_z

    df = pd.read_csv(args.csv).sort_values(KEY + ["step"])
    df["r_mm"] = np.hypot(df["x_mm"], df["y_mm"])
    df["doi"] = df["r_mm"] - Ri
    dep = df[df["process"] != "escape"].copy()        # energy-depositing interactions

    # ---- per-gamma quantities -------------------------------------------------
    Etot = df.groupby(KEY)["e_dep_keV"].sum()
    nblk = dep.drop_duplicates(KEY + ["iz", "iphi"]).groupby(KEY).size()  # distinct blocks hit
    pg = pd.DataFrame({"Etot": Etot})
    pg["nblk"] = nblk.reindex(pg.index).fillna(0).astype(int)             # 0 = pass-through
    pg["contained"] = (pg["Etot"] >= FULL_KEV) & (pg["nblk"] == 1)

    # first / second interaction (the impact and the next deposit)
    dep["rank"] = dep.groupby(KEY).cumcount()
    first = dep[dep["rank"] == 0].set_index(KEY)
    second = dep[dep["rank"] == 1].set_index(KEY)
    first["Etot"] = pg["Etot"]

    # crystal-local face coordinates of the impacts
    phi = np.arctan2(first["y_mm"], first["x_mm"])
    phi_c = (first["iphi"] + 0.5) * dphi
    x_loc = first["r_mm"] * ((phi - phi_c + np.pi) % (2 * np.pi) - np.pi)   # arc offset
    z_c = -H + (first["iz"] + 0.5) * dz
    y_loc = first["z_mm"] - z_c

    # ---- per-event coincidence categories (events with both photons present) --
    cont = pg["contained"].reset_index()
    present = cont.groupby("event_number")["gamma"].nunique()
    both = present[present == 2].index
    ncont = cont[cont["event_number"].isin(both)].groupby("event_number")["contained"].sum()
    n_both, n_cc = len(both), int((ncont == 2).sum())
    cc_events = ncont[ncont == 2].index

    # ---- figure ---------------------------------------------------------------
    fig = plt.figure(figsize=(15, 13))
    fig.suptitle(f"Back-to-back 511 keV — {df['event_number'].nunique()} events with a hit "
                 f"({os.path.basename(args.csv)})", fontsize=13)

    def ax(i):
        return fig.add_subplot(3, 3, i)

    # 1. Edep per gamma (log y: the pass-through spike at 0 dwarfs the photopeak)
    a = ax(1)
    a.hist(pg["Etot"], bins=80, range=(0, 520), color="C0")
    a.axvline(511, color="gray", ls="--", lw=1)
    a.set_yscale("log")
    a.set_xlabel("Edep per gamma [keV]"); a.set_ylabel("gammas")
    a.set_title("Energy deposited per gamma")

    # 2. 3-D impacts on the ring (subsampled)
    a = fig.add_subplot(3, 3, 2, projection="3d")
    s = first.sample(min(3000, len(first)), random_state=0)
    p3 = a.scatter(s["x_mm"], s["y_mm"], s["z_mm"], c=s["Etot"], cmap="viridis", s=3)
    fig.colorbar(p3, ax=a, shrink=0.6, label="Edep [keV]")
    a.set_xlabel("x [mm]"); a.set_ylabel("y [mm]"); a.set_zlabel("z [mm]")
    a.set_title("Impact points on the ring")

    # 3. DOI of the 1st interaction
    a = ax(3)
    a.hist(first["doi"], bins=60, range=(0, wall), color="C3")
    a.set_xlabel("DOI 1st interaction [mm]"); a.set_ylabel("gammas")
    a.set_title("Depth of 1st interaction")

    # 4. DOI of the 2nd interaction
    a = ax(4)
    a.hist(second["doi"], bins=60, range=(0, wall), color="C4")
    a.set_xlabel("DOI 2nd interaction [mm]"); a.set_ylabel("gammas")
    a.set_title("Depth of 2nd interaction")

    # 5. both-gamma-contained fraction (events with both photons hitting)
    a = ax(5)
    cats = [(ncont == 2).sum(), (ncont == 1).sum(), (ncont == 0).sum()]
    bars = a.bar(["both", "one", "neither"], cats, color=["C2", "C7", "C1"])
    for r, c in zip(bars, cats):
        a.text(r.get_x() + r.get_width() / 2, c, f"{c/max(n_both,1):.1%}",
               ha="center", va="bottom", fontsize=9)
    a.set_ylabel("events (both photons hit)")
    a.set_title(f"Contained per coincidence ({n_cc}/{n_both} clean)")

    # 6. Edep of contained gammas
    a = ax(6)
    a.hist(pg.loc[pg["contained"], "Etot"], bins=60, range=(490, 520), color="C2")
    a.axvline(511, color="gray", ls="--", lw=1)
    a.set_xlabel("Edep [keV]"); a.set_ylabel("contained gammas")
    a.set_title("Energy of contained gammas")

    # 7. impact x-y within a crystal (all crystals folded into one face)
    a = ax(7)
    lim = max(args.r_inner_mm * dphi, dz) / 2 * 1.05
    h = a.hist2d(x_loc, y_loc, bins=60, range=[[-lim, lim], [-lim, lim]], cmap="viridis")
    fig.colorbar(h[3], ax=a, label="impacts")
    a.set_xlabel("arc x [mm]"); a.set_ylabel("axial y [mm]")
    a.set_title("Impact within the crystal face"); a.set_aspect("equal")

    # 8. blocks touched per gamma (overspill multiplicity)
    a = ax(8)
    nb = pg.loc[pg["nblk"] >= 1, "nblk"]
    cnt = [(nb == 1).sum(), (nb == 2).sum(), (nb == 3).sum(), (nb > 3).sum()]
    bb = a.bar(["1", "2", "3", ">3"], cnt, color="C5")
    for r, c in zip(bb, cnt):
        a.text(r.get_x() + r.get_width() / 2, c, f"{c/max(len(nb),1):.0%}",
               ha="center", va="bottom", fontsize=9)
    a.set_xlabel("blocks touched"); a.set_ylabel("gammas (>=1 deposit)")
    a.set_title("Overspill: blocks per gamma")

    # 9. ring hit map: φ vs z of impacts (unrolled)
    a = ax(9)
    phi_deg = np.degrees(np.arctan2(first["y_mm"], first["x_mm"]) % (2 * np.pi))
    h = a.hist2d(phi_deg, first["z_mm"], bins=[72, 60],
                 range=[[0, 360], [-H, H]], cmap="viridis")
    fig.colorbar(h[3], ax=a, label="impacts")
    a.set_xlabel("φ [deg]"); a.set_ylabel("z [mm]")
    a.set_title("Ring hit map (unrolled)")

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=120)
    print(f"wrote {args.out}")
    print(f"  gammas with a hit = {len(pg)};  contained = {pg['contained'].mean():.3f}")
    print(f"  both-photon events = {n_both};  clean coincidences (both contained) = "
          f"{n_cc} ({n_cc/max(n_both,1):.3f})")


if __name__ == "__main__":
    main()
