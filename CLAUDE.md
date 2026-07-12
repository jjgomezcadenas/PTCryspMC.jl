# CLAUDE.md — PTCryspMC.jl: PET detector simulation

Orientation for any Claude Code session on this repo.

## Working style

Ask questions plainly, in prose — not as multiple-choice "shopping list" menus (avoid the
AskUserQuestion option-list format). State what you need to know directly.

**Terminal output: never emit blue-rendering text; use bold.** JJ reads responses in a terminal
where markdown links, `inline code` (backticks), and file:line references all render blue. Do not
use any of them in responses — no backticks, no code fences for listings, no clickable links. Emphasise
with **bold** and write paths, filenames, config keys, and code identifiers as plain text.

**Figures: always from a tracked tool.** Every figure that goes into a doc (the LaTeX notes,
`docs/`) must be produced by a checked-in script (a `py/fig_*.py`, or a `scripts/` generator) that
reads the data and writes the image — never with a throwaway inline/ad-hoc command. The script and
the generated image are both committed, so any figure regenerates from the repo. E.g.
`py/fig_activity_profiles.py` → `latex/figures/activity_{truth,vs_tstart}.png`.

## Purpose

Simulate how a PET scanner detects the positron activity left by a proton field, and
write the list of coincidences it would record — the list-mode LOR file (`lors_det.h5`) for a
given detector. The proton transport is done upstream by `ptcryspg4`; this repo begins from the
annihilation points. It runs in two modes, both built & validated (Clinical — a tracer distribution
at a known activity; Proton Activity (API) — a `ptcryspg4` scenario — see "The guide" below).

The simulation runs once per detector, and every detector reads the identical scenario, so
differences come from the detector alone. The analysis that consumes the coincidence list
— range precision, comparing detectors — is separate from the simulation and lives in its own repo
(`CryspBrainSim`, which vendors this repo's `docs/`+`dev/` contracts). It is not described here.

## The guide

Two LaTeX manuals in `docs/` describe the simulator. `PTCryspMC_phys.tex` is the **engine** —
the physics, geometry, transport, and detector response (shared by both modes).
`PTCryspMC_app.tex` is the **application** — the two modes the engine is driven by:
**Clinical** (a radiotracer distribution at a known activity — regions, isotope, acquisition
window) and **Proton Activity (API)** (the positron activity a proton dose leaves, via a
`ptcryspg4` scenario). Both are built & validated.
Read these for the method; this file records the decisions, parameters, and build notes.

## Input — a scenario from ptcrysp-scenarios

The **Proton Activity (API) mode** (`[source].mode="api"`, built & validated) reads a frozen
`ptcryspg4` scenario as its source; the **Clinical mode** (`runs/*.toml` with an activity-driven
tracer distribution) reads none. Point the API driver at one scenario directory by
`[source].scenario_dir` (scenarios are append-only, so the name is the version). A scenario carries:

- `emitters.csv` — annihilation points + isotope per emitter (the spatial source).
- `run_meta.csv` — normalization (target dose, Np per Gy, …).
- `sampling_budget_<scenario>.csv` — measured decays N_j per isotope.
- `isotopes.csv`, `SCHEMA.md` — isotope codes and column meanings.

Stamp the scenario name into every output.

**Source positioning (patient placement).** The frozen scenario's coordinates put the tumour off
the ring centre (native frame: tumour centre z = −25 mm, activity distal edge ≈ −16 mm, activity
peak ≈ −55 mm; the β⁺ activity peaks proximal of the dose). Clinically the patient is positioned so
the target is at isocentre, so `[source].center_on = "distal_edge"` (in `src/scenario.jl`
`load_scenario`) rigidly shifts the emitter pool **and** the phantom in z so the activity distal
edge (the range endpoint, R50) sits at z = 0 — maximising the acceptance of the edge-region LORs.
The scanner stays at the origin; emitters and phantom move together (offset ≈ +16.35 mm, stamped
`source_z_offset_mm`). Default (knob absent) = native frame. **All published survey masters were
produced in the NATIVE (off-centre) frame** — see `dev/status.md` for the pending centred re-run.

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
  by FWHM(E) = a·√(511 keV / E); the time is TOF + the scintillation first-photon arrival
  (jitter = −ln u/(N_det·r0), the photostatistics floor) **plus a per-crystal Gaussian readout
  floor σ_t** (trigger/threshold + SPTR + optics) calibrated to measured CTR — the first-photon
  term alone is far too good (CsI: 49 ps FWHM vs 1.84 ns measured; the full story, calibration
  and references in `latex/scanner_prods.tex`). Keep the per-interaction truth in the singles for
  later CNN work.
- **Selection.** A lower cut on the (smeared) energy: `emin` for the spectrum studies (low, so
  the Compton shoulder shows), raised to the photopeak (`reco_emin_keV` = 511 − 3σ_E per
  crystal) for the reconstruction — applied only in the reco merge, so `lors_truth`/`randoms`
  stay uncut. A symmetric window about 511 keV (half-width scaling with the resolution → a
  sharper detector rejects more scatter) is also available, off by default. A clean coincidence =
  each gamma contained in one block (no overspill); the two hits emerge roughly opposite (observed,
  not enforced). Plus the coincidence-window cut |t1 − t2| ≤ τ, with τ **per scanner** at ~3σ of
  the raw t1−t2 from the CTR study (`scripts/studies/ctr_study.jl`): CsI 1.5 ns, BGO 5 ns.
- **Randoms.** A separate pass over the singles, no re-transport: time-tag each single from the
  activity model (the config's single isotope — F-18 default — in clinic mode; the per-isotope
  scenario activity in API mode), restore the absolute clock, sort, pair cross-annihilation singles within the
  coincidence window τ (a few ns). Rate ≈ 2τS², front-loaded, at most a few percent of the trues.
- **Run once.** Transport runs once → trues + singles. Randoms, and any re-realization the
  downstream analysis needs, are cheap operations on the singles list.

## Current detectors (scope)

Monolithic crystals with continuous 3-D readout, three systems in `data/materials.json`
(each entry carries yield, decay mixture, eres_a, pde, σ_t, σ_xyz — the crystal name pulls
everything; configs override only deliberately):

| | CsI (cryo, ~100 K plateau) | BGO_195K (dry ice) | BGO_77K (LN2) |
|---|---|---|---|
| yield (ph/MeV) / decay (ns) | 100k / 800 | 14k / 1768 | 29k / 1400+8700 (26/74%) |
| eres FWHM @511 → reco cut | 6% → 472 keV | 15% → 413 keV | 10% → 446 keV |
| σ_t per gamma / window τ | 0.35 ns / 1.5 ns | 1.1 ns / 5 ns | 1.1 ns / 5 ns |

Shared: pde 0.40 (Hamamatsu S14160-6050HS), σ_xyz 3.5 mm FWHM per axis (= 1.486 mm σ; the
CRYSP SiPM measurement, made on 2X₀ crystals — hence the depth standard below). **Depth = 2X₀**
(CsI 3.72 cm, BGO 2.236 cm; same transverse block everywhere; 3X₀ variants kept for a future
thick-crystal study). Only BGO_195K + CsI are produced — the two BGO temperatures are
timing-degenerate (yield loss ↔ faster decay), so 195K is the cheaper representative; findings
+ calibration in `latex/scanner_prods.tex`. CsI(Tl) (7%) is still to add; pixelated detectors
(LYSO, standard BGO) report a fixed crystal rather than a continuous position, and come later.

**Scanner-geometry survey.** Beyond the 1 m flagship, a family of closed rings varies radius
(200–387 mm) and AFOV (30–102 cm) at fixed block/constants/source, to map σ_R vs geometry —
scanners `crysp_<tag>_{bgo,csi}_2x0` (e.g. `crysp_chs` = compact head scanner, R200 → ~30 cm
bore; `crysp_r35_50cm` = R350 × 50 cm), each in both crystals. The block face is fixed, so the
φ×z counts quantize (radius, AFOV). Radius beats length for a head target (the CHS outperforms
every R≥300 ring below 1 m at a fraction of the crystal); purity is geometry-independent. Full
survey + acceptance table + the products layout: `latex/scanner_prods.tex`.

## Output — `lors_det.h5`

Per accepted coincidence (one LOR): two 3-D hit positions (mm), two energies (keV), two times
(ns), per-gamma phantom-scatter counts (`nscat1`, `nscat2`: 0 clean, 1 single, ≥2 multiple — so a
LOR separates single from multiple scatter), the timing residual, the absolute decay time
`t_decay_s` (s from acquisition start; gamma 1's decay for randoms — lets downstream emulate a
delayed acquisition start as the cut `t_decay ≥ t_start`), and a truth flag (true / scatter /
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
than float CSV, typed/partial fast reads). `lors_truth.h5` is the clean input to the CTR study
(`scripts/studies/ctr_study.jl`, photopeak-conditioned; the older `examine_dt.jl` + `py/plot_dt.py`
remain for dT dumps) that picks the coincidence window τ (set in `[timing].tau_ns`) **per scanner**
at ~3σ of the raw t1−t2: CsI 1.5 ns, BGO 5 ns (the readout floor σ_t dominates, so τ is
crystal-DEPENDENT — the old crystal-independent τ=3 ns predates σ_t; see `latex/scanner_prods.tex`).
The singles carry `t_rel` (TOF + scintillation jitter + Gaussian σ_t, decay-relative, stamped ONCE
in `simulate_source_mt`); both LOR builders reuse it (no jitter recomputed). Then:
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
and the randoms pass — is built and validated in **both** modes. See `dev/status.md` for the current
state and what remains.

**API (Proton Activity) source — BUILT & VALIDATED** (`[source].mode="api"`). Reads a frozen
`ptcryspg4` scenario (emitters + per-isotope decay budget + phantom) and drives the engine from it:
`src/scenario.jl` (reader + `APISource` + `materialize_api_source`), `Ellipsoid` solid,
`G4_BRAIN_ICRP` material, the isotope singles column, per-isotope randoms timing. The historical
8-step build record is `dev/api_plan.md`. Runs export to the `PtCryspProds/` products tree via
`scripts/run/publish_prod.jl` (+ `run_shards.sh`); layout contract `dev/PRODUCTS.md`, downstream
recipe `dev/data_generation_strategy.md`. Multi-region (non-uniform, layered soft-tissue/brain/bone)
head phantoms are deferred and fully scoped in `dev/multiregion_phantom_plan.md` — a self-contained
brief (data, prerequisites, navigator design, tests) for a future instance; the single-region
uniform phantom path is built.

## Tech stack

- **Julia** — the photon transport and geometry, and the coincidence selection. The
  selection runs as a *streaming* pass over the event-ordered stack
  (`scripts/dev/build_coincidences.jl`, and `scripts/build_true_coincidences_from_singles.jl` for
  the singles stack): O(1) memory, no whole-file load, so it scales to the
  large stacks a full run produces (this is why it is Julia, not Python — see the note
  below).
- **Python** — the control plots and the documentation figures (`py/fig_*.py`), and lighter
  downstream analysis. (Scenario reading is in Julia — `src/scenario.jl`.)
- **CSV and HDF5** in and out (CSV for dev/inspection; HDF5 — quantized Int16, compressed —
  for the production singles stack: ~6× smaller, typed/partial fast reads for the write-once,
  read-many singles list).

Originally the selection/coincidence step was slated for Python; it was moved to a Julia
streaming reader for the memory/throughput on large stacks. Hit formation, the energy window,
the randoms pass, and the scenario reader all followed the same Julia approach; Python stays
for plotting and lighter downstream analysis.

## Reference material

- CRYSP detector numbers: Soleti et al. 2024 (arXiv:2406.13598 — cryo CsI constants + the
  measured 1.84 ns CTR that anchors the σ_t calibration; see `latex/scanner_prods.tex`).
- Upstream source method: `ptcryspg4/docs/simulate_pt_pet.tex`, `handoff.tex`.

## Status / next

**`dev/status.md` is the concise current snapshot + the deferred-work / known-nits register —
read it first for "where are we."** In brief:

**Built and validated:** the whole list-mode chain in **both** source modes — the Geant4-style
geometry, the materials / XCOM cross sections, the photon core (`propagate_photon`); source
injection (Clinical activity-driven + API scenario-driven), the block/wheel ring, transport →
singles → `lors_truth` + `randoms` → `lors_det`; hit formation + smear + energy window + two-block
selection; the randoms pass; the `t_decay_s` column (delayed-start emulation downstream); the
three crystal systems with the CTR-calibrated timing (σ_t) and per-crystal σ_xyz/eres/τ; and the
`PtCryspProds/` products export (per scenario/scanner/crystal, + the detector-independent
`truth/` bundle; the ten-shard `crysp_ring_1m/bgo/fast_1Gy` master is a **frozen reference** —
its config uses the removed hybrid "BGO" entry). `Pkg.test` passes (**1024**).

**Produced** (in `PtCryspProds/`, local): the scanner-geometry survey — ten full ten-shard
masters (5 geometries × BGO_195K/CsI) plus 3 geometry pairs at shard 0; acceptance table and
findings in `dev/status.md` + `latex/scanner_prods.tex`. **Remaining** (mechanical or
downstream): promoting the remaining survey singles to masters on demand; the Documenter
doc-site; and the deferred-but-scoped engine gates (multi-region phantom, open dual-head, mixed
crystal, standalone optical MC). See `dev/status.md` for the full register.

Downstream (separate repo, `CryspBrainSim`): reconstruction, range precision, detector comparison.
