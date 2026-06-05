# CLAUDE.md — PTCryspMC.jl: PET detector simulation and comparison

Orientation for any Claude Code session on this repo.

## Purpose

Simulate candidate PET detectors and rank them for **in-room proton-therapy range
verification**. The figure of merit is **σ(range)** — the precision of the
recovered distal fall-off — at realistic, photon-starved statistics. The reference
detector is the CRYSP family (cryogenic CsI), compared against LYSO and BGO, but
the code is not tied to one design.

The positron-emitter **source is produced upstream** by the `ptcryspg4` repo
(Geant4 proton transport + time-decay handoff) and frozen in the
`ptcrysp-scenarios` data repo. This repo reads a scenario and runs the per-detector
study. It never runs proton transport.

## Input — a scenario from ptcrysp-scenarios

Clone `ptcrysp-scenarios` and point the code at one scenario directory by name
(config or env var; scenarios are append-only, so the name is the version). A
scenario carries:

- `emitters.csv` — annihilation points (the source) + isotope per emitter.
- `run_meta.csv` — normalization (target dose, Np per Gy, …).
- `sampling_budget_<scenario>.csv` — measured decays N_j per isotope (for 1 Gy at a
  timing scenario).
- `isotopes.csv`, `SCHEMA.md` — isotope codes and column meanings.

Read the scenario's own `SCHEMA.md`. Stamp the scenario name into every output.

## Architecture

```
scenario (emitters + budget)
  [B] Analytic detector MC   annihilation events → detector response → coincidences
          ⟹  coincidences_<config>.csv
  [C] Reconstruction         coincidence list → image → σ(range)        DEFERRED
```

Runs **once per detector**; every detector consumes the identical source.

## The analytic detector (Stage B)

Not Geant4 — an analytic γ-transport Monte Carlo (the LXeMC core, vendored and
stripped of its 0νββ apparatus). The geometry is simple (a phantom cylinder + a
crystal ring) and the 511 keV physics is easy to do directly: faster, transparent,
no per-detector Geant4 build.

- **Source injection.** Sample annihilation points from `emitters.csv` (with
  replacement) to the per-isotope budget; emit two back-to-back 511 keV γ,
  isotropic, with acollinearity (~0.5° FWHM). Positron range is already in the ANH
  point.
- **Phantom (treated as water).** Track each γ by Klein–Nishina Compton +
  photoelectric. Compton **must be tracked** (degraded energy + new direction,
  multiple scatter): the patient-scatter background is what the energy window
  discriminates against.
- **Crystal (the ring).** A cylindrical shell (inner R, thickness, axial length L,
  material — all parameters). The interaction depth gives the **DOI**; tracking
  Compton in the crystal gives the photopeak-vs-escape split (differs by material).
  Ring acceptance is a ray–cylinder intersection.
- **Singles → coincidences.** The transport emits a **singles** list (per-module
  position, energy, isotope, scatter flag). A sorter then owns time and number:
  assign each single a time from the activity curve, pair within a coincidence
  window τ (5 ns; randoms ~1%), apply the energy window, label true/scatter/random,
  handle multiples and randoms. Smearing: σ_E (6.3% FWHM at 511 keV, scaling as
  √(511/E)), position/DOI (1.7 mm), σ_t. This split keeps the transport spatial and
  puts all timing/number in the sorter (Python).
- **Backgrounds.** Randoms (from the time structure) and patient scatter; for LYSO,
  a ¹⁷⁶Lu volumetric singles source in the crystal.

## Detector matrix

**CRYSP baseline** (Soleti 2024 / `crysp_for_ht.tex`): ring Ø 77.4 cm, AFOV
102.4 cm, monolithic crystals 48×48×37 mm, **6.3 % FWHM** energy resolution
@ 511 keV, **1.7 mm** 3-D resolution incl. DOI, NEMA sensitivity ~120 kcps/MBq,
~1 µs decay, **no TOF**, **no intrinsic radioactivity**.

**Comparators / variants:** LYSO (~10 % E-res, **¹⁷⁶Lu** intrinsic background —
inject as volumetric source), BGO (poor E-res, no intrinsic activity, high
density); CRYSP variants: room-temperature CsI(Tl), cryogenic BGO,
BGO-core/CsI-wing hybrid.

**Discriminators to sweep:** AFOV/solid-angle (acceptance ≈ L/√(L²+R²)), energy
resolution (scatter rejection on oblique LORs), DOI (parallax), intrinsic
radioactivity (LYSO only), TOF (treated as low-value — show diminishing returns).

## Figure of merit — σ(range)

The budget gives the expected measured decays N_j (per isotope, for a dose and
timing). σ(range) is the spread of the recovered range over **Z Poisson
realizations** of that budget: draw M_j ~ Poisson(N_j), form the source, run the
detector + reconstruction, fit the distal fall-off → R_z; σ(range) = std{R_z} over
Z ~ 100 realizations. For tractability the detector is simulated **once** to a
master coincidence list; the Z realizations are Poisson resamples of that list and
only reconstruction repeats. The realization draws are the old `budget_gen.py`
logic from `ptcryspg4`, which lives here.

## Reconstruction (Stage C) — deferred

List-mode, DOI-aware, optionally TOF. Custom list-mode MLEM/OSEM, or CASToR, or
STIR. Consumes only the coincidence list — never the detector truth.

## Output — `coincidences_<config>.csv`

One row per accepted coincidence: two DOI-resolved 3-D hits (mm), two energies
(keV), two timestamps (ns), and an optional truth flag (0 true, 1 scatter,
2 random). Companion `_meta`: detector config, geometry, energy/time windows, the
scenario read, the budget, realization index, seed.

## Units

Positions **mm**, energies **keV**, times **ns**, dose **Gy** — same as the
upstream files.

## Tech stack

- **Julia** — the transport + geometry (vendored LXeMC core: geometry, NIST cross
  sections, samplers, photon transport).
- **Python** — orchestration: realizations, the coincidence sorter, σ(range),
  plots (mirrors LXeMC's Julia-transport + Python-analysis split).
- **CSV** in and out. Reads scenarios from `ptcrysp-scenarios`.

## Reference material

- `crysp_for_ht.tex` / Soleti et al. 2024 — CRYSP design and detector numbers.
- LXeMC (`~/Projects/XeMC/lxe_mc/LXeMC/design/`): `lxemc.tex`,
  `tracking_and_transport.md`, `geometry_v3.md` — the transport core being vendored.
- CASToR (Merlin 2018) / STIR — reconstruction options.
- The upstream source method: `ptcryspg4/docs/handoff.tex` and
  `ptcryspg4/docs/simulate_pt_pet.tex`.

## Status / next

Not built yet. The source side (`ptcryspg4`) is done and the standard scenario
`head_sobp_1e7` is frozen in `ptcrysp-scenarios`. Build order:

1. Read a scenario; inject back-to-back 511 keV γ from sampled ANH points.
2. Vendor the LXeMC transport core; phantom + crystal-ring geometry (parameters).
3. Transport → singles; the coincidence sorter (time, τ, windows, truth, randoms).
4. Detector configs (CRYSP / LYSO + ¹⁷⁶Lu / BGO) and the discriminator sweeps.
5. σ(range): Z realizations, resample, fit the distal fall-off.
6. Reconstruction (Stage C).
