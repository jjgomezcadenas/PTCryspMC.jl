# CLAUDE.md — PTCryspMC.jl: PET detector simulation

Orientation for any Claude Code session on this repo.

## Working style

Ask questions plainly, in prose — not as multiple-choice "shopping list" menus (avoid the
AskUserQuestion option-list format). State what you need to know directly.

## Purpose

Simulate how a PET scanner detects the positron activity left by a proton field, and
write the list of coincidences it would record — the list-mode LOR file (`lors_det.h5`) for a
given detector. The proton transport is done upstream by `ptcryspg4`; this repo begins from the
annihilation points. It runs in two modes (Clinical — a tracer distribution at a known activity,
built & validated; Proton Activity — a `ptcryspg4` scenario, the target — see "The guide" below).

The simulation runs once per detector, and every detector reads the identical scenario, so
differences come from the detector alone. The analysis that consumes the coincidence list
— range precision, comparing detectors — is separate from the simulation, deferred, and
may live in its own repo. It is not described here.

## The guide

Two LaTeX manuals in `docs/` describe the simulator. `PTCryspMC_phys.tex` is the **engine** —
the physics, geometry, transport, and detector response (shared by both modes).
`PTCryspMC_app.tex` is the **application** — the two modes the engine is driven by:
**Clinical** (a radiotracer distribution at a known activity — regions, isotope, acquisition
window: the validated, currently running mode) and **Proton Activity (API)** (the positron
activity a proton dose leaves, via a `ptcryspg4` scenario: the eventual target, not yet built).
Read these for the method; this file records the decisions, parameters, and build notes.

## Input — a scenario from ptcrysp-scenarios

This is the **Proton Activity (API) mode** — the eventual target, not yet implemented. What
currently runs and is validated is the **Clinical mode** (`runs/*.toml` with an activity-driven
tracer distribution; see `docs/PTCryspMC_app.tex`). When built, the scenario will be the source
for the API mode:
clone `ptcrysp-scenarios` and point the code at one scenario directory by name (config or env var;
scenarios are append-only, so the name is the version). A scenario carries:

- `emitters.csv` — annihilation points + isotope per emitter (the spatial source).
- `run_meta.csv` — normalization (target dose, Np per Gy, …).
- `sampling_budget_<scenario>.csv` — measured decays N_j per isotope.
- `isotopes.csv`, `SCHEMA.md` — isotope codes and column meanings.

Stamp the scenario name into every output.

## How the simulation works (decisions)

A fast, photon-only Monte Carlo. Full method in `docs/PTCryspMC_phys.tex`; the decisions:

- **Input photons.** Draw N_j annihilation points per isotope from the scenario (with
  replacement); each emits two back-to-back 511 keV photons with ~0.5° FWHM acollinearity.
  Order of 10⁸ decays for a 1 Gy field. The positron range is already in the annihilation
  point (set upstream).
- **One annihilation at a time.** Pileup is negligible — the rate spreads over ~10³ crystal
  modules, ~10⁻⁴ occupancy even with a slow crystal — so annihilations are simulated
  independently. Each yields, if both photons are detected, one same-annihilation
  coincidence (true/scatter); every detected hit also goes to a singles list.
- **Transport.** Follow only the photons (electrons deposit locally; no scintillation). Free
  path from XCOM cross sections; Compton (Klein–Nishina) + photoelectric; pair closed at
  511 keV. Phantom treated as water (Compton scatter = the patient-scatter background);
  crystal ring as a cylindrical shell.
- **Geometry: structured block/wheel grid.** The ring is segmented in φ (blocks) and z
  (wheels) as a structured grid on the shell — block index by arithmetic, boundary distances
  closed-form — so the cost is independent of block count. Dead gaps skipped for now.
  Overspill (scatter into a neighbour block) and energy non-containment (scattered photon
  escapes) come out of the transport.
- **Hits.** Per block: energy = sum, position = **first interaction point** (not the
  centroid — it is the LOR point and the target a CNN recovers when separating 1- from
  2-Compton events), time = earliest. Smear the position by σ_xyz (incl. DOI) and the energy
  by FWHM(E) = a·√(511 keV / E); the time is the scintillation first-photon arrival (TOF +
  jitter = −ln u/(N_det·r0)), a physical distribution rather than a Gaussian σ_t. Keep the
  per-interaction truth in the singles for later CNN work.
- **Selection.** A lower cut on the (smeared) energy: `emin` for the spectrum studies (low, so
  the Compton shoulder shows), raised to the photopeak region (`reco_emin_keV`) for the
  reconstruction. A symmetric window about 511 keV (half-width scaling with the resolution → a
  sharper detector rejects more scatter) is also available, off by default. A clean coincidence =
  each gamma contained in one block (no overspill); the two hits emerge roughly opposite (observed,
  not enforced). Plus the coincidence-window cut |t1 − t2| ≤ τ.
- **Randoms.** A separate pass over the singles, no re-transport: time-tag each single from the
  activity model (the config's single isotope — F-18 default — in clinic mode; the per-isotope
  scenario activity in API mode), restore the absolute clock, sort, pair cross-annihilation singles within the
  coincidence window τ (a few ns). Rate ≈ 2τS², front-loaded, at most a few percent of the trues.
- **Run once.** Transport runs once → trues + singles. Randoms, and any re-realization the
  downstream analysis needs, are cheap operations on the singles list.

## Current detectors (scope)

Monolithic crystals with continuous 3-D readout: pure CsI (a = 5% FWHM @ 511 keV) and
cryogenic BGO (10%), both validated; position/DOI σ_xyz ~1.7 mm (CRYSP). CsI(Tl) (7%) is
still to add. Pixelated detectors (LYSO, standard BGO, a = 15–20%) report a fixed crystal
rather than a continuous position, and come later.

## Output — `lors_det.h5`

Per accepted coincidence (one LOR): two 3-D hit positions (mm), two energies (keV), two times
(ns), per-gamma phantom-scatter counts (`nscat1`, `nscat2`: 0 clean, 1 single, ≥2 multiple — so a
LOR separates single from multiple scatter), the timing residual, and a truth flag (true / scatter /
random). The HDF5 root attributes carry the provenance: detector config, geometry, energy and time
windows, the run tag, and the seed.

The full column schema for every output file (singles + LORs, column · type · unit · meaning, the
quantization, the truth code, the provenance attrs) is **`docs/SCHEMA.md`**, generated from the code
by `scripts/gen_schema.jl` (column names/types introspected from `singles_columns`/`coinc_columns`,
units/meaning from the co-located `singles_doc`/`coinc_doc` maps). A test fails if it drifts — so
regenerate it after any schema change rather than hand-editing.

## Code layout (`src/`)

The photon transport and geometry are photon-only:

- `geometry.jl` — a Geant4-style hierarchy: `Solid` (`Cylinder`, `Box`, `Sphere`, `CylShell`),
  `LogicalVolume` (solid + material), `PhysicalVolume` (a placed logical volume), `Scanner` (the
  ring + its block/wheel grid), with `is_inside`, `distance_to_entry`, `distance_to_exit`.
  `load_geometry` reads the world from `geometry/geometry.json` (named sections: `world`,
  `phantom`, `scanner`) into a `Geometry` container; `load_solid` is the shape factory.
- `nist_data.jl`, `materials.jl` — the XCOM loader and `sigma_macro(material, E)` (the
  macroscopic Compton/photoelectric/pair cross sections); `load_material` /
  `load_materials`. The water and crystal tables (`data/xcom_{water,CSI,BGO}.csv`) are in place,
  with the per-crystal scintillation properties (light yield, decay, eres, PDE); tissue and LYSO
  still to add.
- `sampling.jl`, `transport.jl` — the Compton/photoelectric samplers and the photon step
  loop `propagate_photon` (through a `PhysicalVolume`); only the photon path is followed,
  with local energy deposition.

Other dirs: `geometry/` (JSON world), `data/` (materials + XCOM), `scripts/` (the **production
chain** — `simulate_source_mt.jl`, `build_true_coincidences_from_singles.jl`,
`build_randoms_from_singles.jl`, `reco_lors.jl`; subdirs `scripts/dev/`
= the full-stack dev chain `simulate_phantom.jl` + `build_coincidences.jl`, `scripts/studies/`
= one-off explorations, `scripts/tests/` = QA/benchmark scripts, `scripts/run/` = launchers
(`run_matrix.sh` fans the dev chain across configs → `output/`; `run_prod.sh` runs the production
chain per NAMED config → `prod/` (explicit-only, HDF5 + pinned nchunks for reproducibility))), `py/` (Python plotters),
`runs/` (TOML run configs, tracked), `test/`. Three
gitignored output trees, one per script category (each holding only directories): `output/`
(dev chain — `<tag>/` run cases + `control_plots/`), `prod/` (production chain — `<tag>/`
singles + LORs), `studies/` (`scripts/studies/` outputs, by topic: `lifetime/`, `xsections/`,
`b2b/`, `crystal/`, `phantom/`).

The **production chain** writes under **`prod/<tag>/`** (separate from the dev `output/<tag>/`;
base from `[output].prod_dir`, default `prod`): `scripts/simulate_source_mt.jl` (multi-threaded,
singles-only via the allocation-free `navigate_single_photons`) → `prod/<tag>/singles.{csv,h5}`,
then `scripts/build_true_coincidences_from_singles.jl` (reads singles either-format, fills `GammaAcc`
directly via the shared `src/coincidences.jl`) → `prod/<tag>/lors_truth.h5` — the same-annihilation
coincidences (true + scatter), **truth-only**: no detector smearing, no energy cut, but carrying the
per-gamma timestamps `t1,t2` and the residual `dt = |Δt0| − TOF_diff`. HDF5 stores quantized Int16
columns (0.1 mm / 0.1 keV — lossless at detector resolution), chunked + shuffle+deflate (~6× smaller
than float CSV, typed/partial fast reads). `lors_truth.h5` is the clean input to the DT study
(`examine_dt.jl` → `studies/dt/<tag>_dt.csv`, `py/plot_dt.py`) that picks the coincidence window τ
(set in `[timing].tau_ns`; τ=3 ns for CsI & BGO — the dT tail is TOF/geometry-set, crystal-independent).
The singles carry `t_rel` (TOF + scintillation jitter, decay-relative, stamped ONCE in
`simulate_source_mt`); both LOR builders reuse it (no jitter recomputed). Then:
`build_randoms_from_singles.jl` → `randoms.h5` (truth=2): restore the absolute clock per single
(`event_time(ev)·1e9 + t_rel`, in memory only), sort, pair cross-decay singles within τ
(`src/randoms.jl pair_randoms`); and `reco_lors.jl` → **`lors_det.h5`** (`mode=det`,
`has_randoms=true`): stream truth ∪ randoms through smear (σ_xyz, eres-from-crystal) + energy-select
+ DT cut (`|t1−t2|≤τ`) + the truth flag (true/scatter/random) — a streaming concat, NOT time-sorted
(order-independent for recon; a chronological stream is recomputable from `(event,t1)` on demand).
`scripts/tests/check_singles.jl` / `check_lors.jl` validate the stacks; `diff_singles.jl` compares two.

The LOR-generation pipeline (`simulate_phantom.jl` → `build_coincidences.jl` →
`plot_coincidences.py`) is **TOML-config driven**: a `runs/<tag>.toml` (sections
`[geometry] [source] [transport] [detector] [timing] [output]`) is the parameter source of truth
and the run's provenance. The `tag` (= config filename, overridable) names the per-run
output dir `output/<tag>/`, which holds the stack, the coincidence file(s), the plot, and a
copy of the config. `src/config.jl` reads it (`read_config`, `run_tag`, `cfg_get`).

The full chain — source injection, the block/wheel grid, transport, hit formation + selection,
and the randoms pass — is built and validated in Clinical mode. See `dev/status.md` for the current
state and what remains (the Proton Activity (API) scenario source; the Documenter doc-site).

**API (Proton Activity) source — in progress on branch `api-scenario`.** Reads a frozen
`ptcryspg4` scenario (emitters + per-isotope decay budget + phantom) and drives the engine from it.
`dev/api_plan.md` is the step-by-step build plan (Ellipsoid solid, brain material, phantom loader,
scenario reader, APISource + Poisson materialization, isotope column, driver). Multi-region
(non-uniform, layered soft-tissue/brain/bone) head phantoms are deferred and fully scoped in
`dev/multiregion_phantom_plan.md` — a self-contained brief (data, prerequisites, navigator design,
tests) for a future instance; the single-region uniform phantom path is built.

## Tech stack

- **Julia** — the photon transport and geometry, and the coincidence selection. The
  selection runs as a *streaming* pass over the event-ordered stack
  (`scripts/dev/build_coincidences.jl`, and `scripts/build_true_coincidences_from_singles.jl` for
  the singles stack): O(1) memory, no whole-file load, so it scales to the
  large stacks a full run produces (this is why it is Julia, not Python — see the note
  below).
- **Python** — scenario reading (API mode, to come), the control plots and the
  documentation figures (`py/fig_*.py`), and lighter downstream analysis.
- **CSV and HDF5** in and out (CSV for dev/inspection; HDF5 — quantized Int16, compressed —
  for the production singles stack: ~6× smaller, typed/partial fast reads for the write-once,
  read-many singles list).

Originally the selection/coincidence step was slated for Python; it was moved to a Julia
streaming reader for the memory/throughput on large stacks. Hit formation, the energy
window and the randoms pass will follow the same Julia-streaming approach; Python stays
for scenario reading and plotting.

## Reference material

- CRYSP detector numbers: Soleti et al. 2024.
- Upstream source method: `ptcryspg4/docs/simulate_pt_pet.tex`, `handoff.tex`.

## Status / next

**`dev/status.md` is the concise current snapshot + the deferred-work / known-nits register —
read it first for "where are we."** In brief:

**Built and validated:** the foundations — the Geant4-style geometry, the materials / XCOM
cross sections, and the photon physics core (`propagate_photon`) — plus the first result:
511 keV photons transported through the water phantom to a photon stack, reproducing
Beer–Lambert (unscattered fraction 0.215 vs 0.216). `Pkg.test` passes.

**Done (1–5):** scenario injection, phantom + block/wheel ring, transport → singles +
same-annihilation coincidences, hit formation + smear + energy window + two-block selection,
**and the randoms pass** — the full list-mode chain `singles → lors_truth + randoms → lors_det`
(per-photon time `t_rel` stamped once in the singles; trues reuse it, randoms restore the absolute
clock to pair cross-decay singles within τ; reco merges + smears + cuts + flags). See the chain
description above.

**Remaining:**

6. The monolithic detector configs (CsI/CsI(Tl)/BGO sweep) + benchmark/MT for the 10⁸ production
   (external sort / τ-bin streaming where the in-memory passes don't scale).

Downstream (separate, deferred): range precision and detector comparison; reconstruction.
