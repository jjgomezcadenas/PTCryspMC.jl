# CLAUDE.md — PTCryspMC.jl: PET detector simulation

Orientation for any Claude Code session on this repo.

## Purpose

Simulate how a PET scanner detects the positron activity left by a proton field, and
write the list of coincidences it would record. The source comes from upstream
(`ptcryspg4`, frozen in `ptcrysp-scenarios`); this repo reads a scenario and produces
`coincidences_<config>.csv` for a given detector. It never runs proton transport.

The simulation runs once per detector, and every detector reads the identical scenario, so
differences come from the detector alone. The analysis that consumes the coincidence list
— range precision, comparing detectors — is separate from the simulation, deferred, and
may live in its own repo. It is not described here.

## The guide

`docs/pet_simulation.tex` explains the simulation: the idea, the input photons, how the
fast Monte Carlo works (interactions, geometry, transport), and what it produces. Read it
first. This file records the decisions, parameters and build notes, and refers to the
guide for the method.

## Input — a scenario from ptcrysp-scenarios

Clone `ptcrysp-scenarios` and point the code at one scenario directory by name (config or
env var; scenarios are append-only, so the name is the version). A scenario carries:

- `emitters.csv` — annihilation points + isotope per emitter (the spatial source).
- `run_meta.csv` — normalization (target dose, Np per Gy, …).
- `sampling_budget_<scenario>.csv` — measured decays N_j per isotope.
- `isotopes.csv`, `SCHEMA.md` — isotope codes and column meanings.

Stamp the scenario name into every output.

## How the simulation works (decisions)

A fast, photon-only Monte Carlo. Full method in `docs/pet_simulation.tex`; the decisions:

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
  2-Compton events), time = earliest. Smear by σ_xyz (incl. DOI), σ_t, and an
  energy-dependent FWHM(E) = a·√(511 keV / E). Keep the per-interaction truth in the singles
  for later CNN work.
- **Selection.** Energy window symmetric about 511 keV, half-width ~2 FWHM (tunable; its
  width scales with the resolution → a sharper detector rejects more scatter). A clean
  coincidence = exactly two roughly opposite blocks touched; overspill into a third block,
  an under-energy escape, or a miss fails it.
- **Randoms.** A separate pass over the singles, no re-transport: time-tag each single from
  its isotope's activity curve, sort, pair cross-annihilation singles within the coincidence
  window τ (a few ns). Rate ≈ 2τS², front-loaded, at most a few percent of the trues.
- **Run once.** Transport runs once → trues + singles. Randoms, and any re-realization the
  downstream analysis needs, are cheap operations on the singles list.

## Current detectors (scope)

Monolithic crystals with continuous 3-D readout: pure CsI (a = 5% FWHM @ 511 keV), CsI(Tl)
(7%), cryogenic BGO (10%); position/DOI σ_xyz ~1.7 mm (CRYSP). Pixelated detectors (LYSO,
standard BGO, a = 15–20%) report a fixed crystal rather than a continuous position, and
come later.

## Output — `coincidences_<config>.csv`

Per accepted coincidence: two 3-D hit positions (mm), two energies (keV), two times (ns),
and a truth flag (true / scatter / random). Companion `_meta`: detector config, geometry,
energy and time windows, scenario name, seed.

## Code layout (`src/`)

The photon transport and geometry are photon-only:

- `geometry.jl` — a Geant4-style hierarchy: `Solid` (`Cylinder`; `CylShell` to come),
  `LogicalVolume` (solid + material), `PhysicalVolume` (a placed logical volume), with
  `is_inside`, `distance_to_entry`, `distance_to_exit`. `load_geometry` reads the world
  from `geometry/geometry.json` (named section per component, `phantom` so far) into a
  `Geometry` container; `load_solid` is the shape factory.
- `nist_data.jl`, `materials.jl` — the XCOM loader and `sigma_macro(material, E)` (the
  macroscopic Compton/photoelectric/pair cross sections); `load_material` /
  `load_materials`. The water table (`data/xcom_water.csv`, 10 keV–10 MeV) is in place;
  tissue and the crystals (CsI, BGO, LYSO) still to add.
- `sampling.jl`, `transport.jl` — the Compton/photoelectric samplers and the photon step
  loop `propagate_photon` (through a `PhysicalVolume`); only the photon path is followed,
  with local energy deposition.

Other dirs: `geometry/` (JSON world), `data/` (materials + XCOM), `scripts/` (the **production
chain** — `simulate_source_mt.jl`, `build_coincidences_from_singles.jl`; subdirs `scripts/dev/`
= the full-stack dev chain `simulate_phantom.jl` + `build_coincidences.jl`, `scripts/studies/`
= one-off explorations, `scripts/tests/` = QA/benchmark scripts, `scripts/run/` = parallel
launchers), `py/` (Python plotters), `runs/` (TOML run configs, tracked), `output/` (per-run
`output/<tag>/` results, gitignored), `test/`.

The **production chain** writes under **`prod/<tag>/`** (separate from the dev `output/<tag>/`;
base from `[output].prod_dir`, default `prod`): `scripts/simulate_source_mt.jl` (multi-threaded,
singles-only via the allocation-free `navigate_single_photons`) → `prod/<tag>/singles.{csv,h5}`,
then `scripts/build_coincidences_from_singles.jl` (reads singles either-format, fills `GammaAcc`
directly via the shared `src/coincidences.jl`) → `prod/<tag>/lors_{truth,det}.h5`. HDF5 stores
quantized Int16 columns (0.1 mm / 0.1 keV — lossless at detector resolution), chunked +
shuffle+deflate (~6× smaller than float CSV, typed/partial fast reads). The LOR list is the
list-mode deliverable (truth ∈ {true,scatter}; random=2 + a time-ordered merge are Step 5,
`has_randoms=false` until then). `scripts/tests/check_singles.jl` / `check_lors.jl` validate the
stacks; `diff_singles.jl` compares two.

The LOR-generation pipeline (`simulate_phantom.jl` → `build_coincidences.jl` →
`plot_coincidences.py`) is **TOML-config driven**: a `runs/<tag>.toml` (sections
`[geometry] [source] [transport] [detector] [output]`) is the parameter source of truth
and the run's provenance. The `tag` (= config filename, overridable) names the per-run
output dir `output/<tag>/`, which holds the stack, the coincidence file(s), the plot, and a
copy of the config. `src/config.jl` reads it (`read_config`, `run_tag`, `cfg_get`).

Still to write: source injection, the block-grid index and plane-crossing distances, hit
formation + selection, the randoms pass.

## Tech stack

- **Julia** — the photon transport and geometry, and the coincidence selection. The
  selection runs as a *streaming* pass over the event-ordered stack
  (`scripts/dev/build_coincidences.jl`, and `scripts/build_coincidences_from_singles.jl` for
  the singles stack): O(1) memory, no whole-file load, so it scales to the
  large stacks a full run produces (this is why it is Julia, not Python — see the note
  below).
- **Python** — read a scenario, the control plots, and lighter downstream analysis.
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

`dev/dev_steps.md` is the running build log (target, progress, next steps); consult it for
the current state. In brief:

**Built and validated:** the foundations — the Geant4-style geometry, the materials / XCOM
cross sections, and the photon physics core (`propagate_photon`) — plus the first result:
511 keV photons transported through the water phantom to a photon stack, reproducing
Beer–Lambert (unscattered fraction 0.215 vs 0.216). `Pkg.test` passes.

**Remaining build order:**

1. Read a scenario; inject the back-to-back 511 keV pairs.
2. Phantom + block/wheel ring geometry (`CylShell` + the structured grid).
3. Transport → singles + same-annihilation coincidences (multi-volume navigator).
4. Hit formation (first interaction, smear, energy window) and the two-block selection.
5. The randoms pass.
6. The monolithic detector configs.

Downstream (separate, deferred): range precision and detector comparison; reconstruction.
