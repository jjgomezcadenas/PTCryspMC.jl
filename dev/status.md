# Project status — PTCryspMC.jl

A one-page snapshot of where the simulation stands, plus the **deferred-work register**.
Companion docs: the *method* is in `docs/PTCryspMC_phys.tex` (engine) and `docs/PTCryspMC_app.tex` (modes); the
decisions + code layout in `CLAUDE.md`.

_Last updated: 2026-07-09._

## What it does

A fast, photon-only Monte Carlo: from a **source** (a clinical tracer distribution, or the positron
activity a proton field leaves) + a detector description, write the list-mode coincidence/LOR file a
PET scanner would record. It never runs proton transport. Both source modes are built & validated.

## The pipeline (built & validated)

```
runs/<tag>.toml
   │
simulate_source_mt.jl  (multi-threaded, alloc-free) ──► prod/<tag>/singles.h5   (+ t_rel, n_scatter)
   ├─ build_true_coincidences_from_singles.jl ──► lors_truth.h5   (true + scatter; t1,t2,dt, nscat1,nscat2)
   ├─ build_randoms_from_singles.jl            ──► randoms.h5      (truth = random)
   └─ reco_lors.jl  (smear + energy + DT cut + flag) ──► lors_det.h5   (the list-mode deliverable)
```

- **Scatter multiplicity.** Each photon carries a phantom-scatter **count** `n_scatter` (Compton
  interactions in the phantom); the LOR carries both as `nscat1`/`nscat2`, so a coincidence separates
  true (`==0`), single (`nscat1+nscat2==1`) and multiple (`≥2`) scatter for scatter correction.

- **Timing model.** Each single carries `t_rel` = TOF + scintillation first-photon jitter
  (`jitter = −ln u / (N_det·r0)`, `N_det = yield·E·pde`, `r0 = Σ wₖ/τₖ`) **+ a per-crystal
  Gaussian readout floor σ_t** (`sigma_t_ns` on the crystal: trigger/threshold + SPTR + optics —
  the first-photon term alone is the photostatistics floor, ~40× too good for CsI vs the
  measured 1.84 ns CTR of Soleti et al. 2024; calibration + references in `latex/ctr_note.tex`),
  stamped **once** at generation, **relative to the decay** so it stays small (Float32).
  Absolute time = `event_time(ev)·1e9 + t_rel`, reconstructed only where randoms need a common clock.
- **Absolute decay time on every LOR.** All three LOR files carry `t_decay_s` (Float32 s, zero =
  acquisition start; gamma 1's decay for randoms, the `x0` convention) — the `event_time(ev)` of
  the annihilation, so downstream can emulate any delayed acquisition start as the pure cut
  `t_decay ≥ t_start` (CryspBrainSim request `upstream_request_lor_decay_time.md`; ~3.5 B/row
  compressed, ~10% file growth). Attr `t_decay_zero = "acquisition_start"`.
- **Coincidence window τ per scanner**, at ~3σ of the raw `t1−t2` from the CTR study
  (`scripts/studies/ctr_study.jl`, photopeak-conditioned — uncut few-keV deposits have
  1/E-divergent jitter): **CsI 1.5 ns** (raw σ 0.51 ns, keeps 99.7% of trues), **BGO 5 ns**
  (raw σ 1.77 ns, keeps 99.3%; both temperatures timing-identical). With σ_t dominant, τ is
  crystal-DEPENDENT — the old crystal-independent τ = 3 ns (and `compare_crystal_timing.jl`'s
  rationale) predates the readout floor. Windows anchored: CsI to the Soleti bench CTR, BGO to
  the ~5 ns window of commercial BGO systems (GE Omni Legend).
- **Reco lower energy cut** `reco_emin_keV = 511 − 3σ_E` per crystal (CsI 472 / BGO_77K 446 /
  BGO_195K 413 keV; the spectrum studies keep `emin_keV = 300` to see the Compton shoulder).
  Applied ONLY in the reco merge — `lors_truth`/`randoms` stay uncut. No upper cut (analysis-time).
- **Validated:** `Pkg.test` **1024**; randoms match the analytic `2τS²` (CsI 1248 vs 1291, ratio 0.97;
  the clinical Vacuum/BGO 10⁸ runs at fixed N over 100× in rate — 100 kHz×1000 s: 52694 vs 52708;
  1 MBq×100 s: 579705 vs 578399; 10 MBq×10 s: 5837886 vs 5839021 — all ratio 1.00, so 2τS² holds across
  three orders of magnitude in activity, with randoms reaching ~10% of trues at the 10 MBq end);
  reco acceptance CsI 8.98% / BGO 23.75% (historical, pre-2X₀/pre-σ_t numbers);
  the corrected-residual medians matched the analytic single-photon jitter in the pre-σ_t model,
  and the CTR study validates the calibrated model (CsI CTR 1.17 ns / BGO 4.1 ns FWHM);
  the clinical N-from-activity matches the analytic to the event.
- **Scale.** `scripts/tests/bench_chain.jl`: the full chain at 10⁸ runs **~3 min serial, ~14 GB
  peak** (`build_randoms`, the only N-scaling stage) — fits 48 GB, no rework needed. (The clinical
  Vacuum/BGO 10⁸ run completed the whole chain in 180 s on 16 threads.)

## Source scenarios (the front end)

The engine + chain above are shared; what drives them is the **source**, in one of two scenarios
selected by `[source].mode`:

- **Clinical** (activity-driven) — *built & validated.* A tracer distribution at a known activity:
  one or more `[[source.region]]` (a named geometry `volume`, or an inline `shape` + dims +
  `position_cm`), each with a concentration (`conc_kBq_per_mL`) or a total (`activity_kBq`). An isotope
  (F-18 default → `T½`, `β⁺`) + an acquisition window `[t0_s,t1_s]` fix the rest:
  `N = β⁺·(Σ cᵢVᵢ/λ)(1−e^{−λT})` annihilations, each drawn from a region ∝ its activity `cᵢVᵢ` and
  timed by the decay curve. Structured phantoms (Derenzo, NEMA-IQ) need no geometry change (the
  inserts share the phantom material → the transport is unchanged). The count-driven uniform phantom
  is the pinned-N special case (back-compat). QA: `scripts/tests/check_clinic_regions.jl`. Configs:
  `sphere_water_f18_csi`, the sphere/cylinder × air/water BGO set, and the NEMA quartet
  (`nema_{air,water}_bgo`, `nema_la_{air,water}_bgo`). `run_prod.sh` derives N for clinic
  configs (no `--nevents`) and runs NAMED configs only.
- **Proton Activity (API)** (count-driven) — *built & validated.* A frozen `ptcryspg4` scenario
  supplies the emitters and the per-isotope decay budget; the source materializes as
  `M_j ~ Poisson(N_j·f_inside)`, seeded by `(master_seed, realization)` independent of the transport
  chunking. Full detail in the "API source + products handoff" section below; the historical 8-step
  build record is `dev/api_plan.md`.

## Detector configs

Three crystal systems in `data/materials.json`, each carrying its full response (yield, decay,
eres_a, pde 0.40, σ_t, σ_xyz — configs override only deliberately):

| | CsI (cryo plateau, ~100 K) | BGO_195K (dry ice) | BGO_77K (LN2) |
|---|---|---|---|
| yield / decay | 100k ph/MeV / 800 ns | 14k / 1768 ns | 29k / 1400+8700 ns (26/74%) |
| eres → reco cut | 6% → 472 keV | 15% → 413 keV | 10% → 446 keV |
| σ_t / τ | 0.35 ns / 1.5 ns | 1.1 ns / 5 ns | 1.1 ns / 5 ns |

σ_xyz = 3.5 mm FWHM per axis for all (CRYSP SiPM measurement on 2X₀ crystals). **Depth
standard: 2X₀** (CsI 3.72 cm / BGO 2.236 cm — resolution measured at 2X₀; 3X₀ would degrade it
and raise cost; 3X₀ geometries kept for a thick-crystal study). Scanners named with the depth:
`crysp_ring_1m_{bgo,csi}_{2x0,3x0}`. The BGO temperatures are timing-degenerate → **the two
representative production scanners are BGO_195K and CsI** (`runs/uniform_headep_{bgo195k,csi}_api.toml`;
`bgo77k` config ready but not scheduled). Rationale + findings: `latex/ctr_note.tex`. The old
hybrid "BGO" material is removed; `runs/uniform_headep_bgo_api.toml` is the frozen reference of
the existing `crysp_ring_1m/bgo` master. CsI(Tl) (7%) and the pixelated detectors (LYSO,
standard BGO) are still to add.

## Documentation

- `docs/PTCryspMC_phys.tex` — the engine (physics, geometry, transport, detector response).
- `docs/PTCryspMC_app.tex` — the application (the Clinical and API source scenarios).
- Docstrings + a Literate/Documenter web doc-site — **planned (the next doc task), not built.**

## Deferred work & known nits

From the multi-agent code review. **None are correctness bugs** — all are robustness/clarity.
Each numbered item has a natural trigger that makes it cheap to fold in *then* rather than now;
the corresponding code site carries a short `Deferred:` marker pointing here.

### Deferred fixes (specific trigger)

| # | Where | What | Trigger |
|---|-------|------|---------|
| 1 | `src/geometry.jl` `Scanner.volume` | the field named `volume` (a `PhysicalVolume`) collides with the generic `volume()` (cm³); rename → `pvol`/`shell`. | A small rename sweep (`geometry.jl` + `simulate_source_mt.jl`); do it when next touching the Scanner API. |
| 2 | `src/nist_data.jl` `load_xcom` | the parser assumes 8 clean numeric columns + digit-first data rows; XCOM **K-edge label rows** can shift columns, and empty input throws obscurely. | Adding the **CsI/BGO/LYSO XCOM tables** (which have in-range K-edges) — harden with the real tables to test against. |

_(The former #1 — the `first_photon_jitter` `-log(1-rand)` guard — is now fixed, folded into the `n_scatter` schema regeneration.)_

### Minor nits (cosmetic / robustness, no trigger)

- `block_index` docstring (geometry): φ ∈ [0,2π) never equals 2π (the clamp is FP-defensive); `z < −H` silently clamps to `iz=0`.
- `_prepare_xcom_energy` (nist_data): the K-edge nudge assumes the below-edge row precedes the above-edge row (true for standard XCOM) — note it when crystals arrive.
- `_leaf_reduce` (navigator): the inner Compton loop has no explicit step backstop; it terminates because energy strictly decreases each Compton to the cut — worth a one-line comment.
- `position_cm` (geometry loader): a JSON `position_cm` with ≠3 entries throws a confusing `TypeAssertionError`; add a `length==3` check.
- HDF5 chunk size: a fixed 2²⁰-row chunk regardless of run size (oversized/zero-padded for tiny runs); the `total==0` empty-run path is fine but undocumented.
- `build_randoms` Float32 edge: `t2 = Float32(t_rel+Δ)`, then reco recomputes `|t1−t2|`; Float32 rounding can nudge a `Δ ≈ τ` random just across the cut. Sub-permille.

### Deferred by scope (not from the review)

- Multi-region clinic spatial draw uses dynamic dispatch per event — worth a glance if a clinic run
  ever goes to 10⁸ with many regions.
- Refine the per-crystal **PDE** (0.45 placeholder for CsI and BGO; should differ by emission colour).
- Threshold / CFD timing (the "first photon" model is the leading-edge idealization).
- **Open dual-head geometry** (the CRYSP-open arm of the range-verification study,
  `docs/range_verification_recipe.md`) — postponed. When taken up, the cheap representation is a
  *partial ring*: a `phi_gaps` angular acceptance on the existing `CylShell` `Scanner` (blocks in
  the gap arcs void), which captures the missing angular coverage without new solids/navigation;
  flat two-panel geometry only if the planar detail itself becomes the question. The limited-angle
  penalty materializes only through reconstruction (MLEM — downstream, not in this repo).
- **Range-study analysis migrates out** when the reconstruction (MLEM) repo is created: `git mv`
  `docs/range_verification_recipe.md` + `py/range_endpoint.py` there. What stays HERE is the
  master-production plan the recipe fixes (per geometry: ten independent 10⁸-decay trues-only
  shards, distinct seeds — 10× the top dose point), since it runs on this repo's chain.
- Pixelated detectors report a fixed crystal, not a continuous position. Placement **rotation** transform (only when a volume needs it).

## API source + products handoff — BUILT & VALIDATED (2026-07-05)

- **API (Proton Activity) source** — the second source branch is DONE (`[source].mode="api"`).
  Reads a frozen `ptcryspg4` scenario (phantom from `phantom_regions.csv`, per-isotope emitter
  pools, decay budget), materializes the source (`M_j ~ Poisson(N_j·f_inside)`, seeded by
  `(master_seed, realization)` independent of the transport chunking), drops escaped positrons,
  transports, and writes a self-describing `lors_det.h5`. Full 8-step build + validation in
  **`dev/api_plan.md`**; engine pieces: `Ellipsoid` solid, `G4_BRAIN_ICRP` material + `xcom_brain.csv`,
  `src/scenario.jl` (reader + `APISource` + `materialize_api_source`), isotope singles column,
  per-isotope randoms timing. QA: `scripts/tests/check_scenario.jl`, `check_api_source.jl`,
  `check_api_validation.jl`. **`Pkg.test` 1024.**
- **Products handoff (`PtCryspProds`)** — `scripts/run/publish_prod.jl` exports a run into the tree
  `<scenario>/<scanner>/<crystal>/<budget>_<dose>/lors_shardNNN.h5` (+ shared `scanner_geometry.json`,
  `phantom/`, `README.md`); `scripts/run/run_shards.sh` generates shards sequentially (chain →
  publish → gates → prune all `.h5`). Layout contract in **`dev/PRODUCTS.md`**, the downstream
  recipe (shards vs realizations, `thin_lm`, σ_R) in **`dev/data_generation_strategy.md`**.
- **Truth bundle (`truth/`)** — `publish_prod` also exports a detector-independent `<scenario>/truth/`
  (shared like `phantom/`): the dose-side truth (`depth_dose.csv` → dose-R80, `sobp_layers`, `run_meta`,
  `sampling_budget_<budget>`) plus the derived `activity_profile_<budget>.csv` (binned true activity(z)
  per isotope → activity-R50, on the `depth_dose` z-frame, scaled `N_expected·f_inside` so it composes
  with the shards). Runnable stand-alone (`publish_prod --truth-only`) to backfill. Requested by the
  downstream repo in `CryspBrainSim/dev/upstream_request_truth_bundle.md`.
- **First master produced:** `~/Projects/PtCryspProds/uniform_headep_sobp_1e8/crysp_ring_1m/bgo/
  fast_1Gy/` — 10 shards (realizations 0–9, ΣM = 8.02e8, all validated) + the `truth/` bundle.
  `runs/uniform_headep_bgo_api.toml`.

## Next

1. **CsI arm + other scanners** — a CsI config (copy the BGO one, `crystal_material="CsI"`) →
   `run_shards … 0 9` for the detector comparison at matched shards; head/children `CylShell`
   variants for the geometry comparison. All mechanical now.
2. **Downstream reconstruction / σ_R repo** — consumes `PtCryspProds` (see
   `dev/data_generation_strategy.md`); builds `thin_lm` (Bernoulli, pooled shards), MLEM, the
   depth profile; `py/range_endpoint.py` has `fit_endpoint`/`sigma_R`.
3. **Deferred, scoped:** multi-region (skull/brain/scalp) head phantom — `dev/multiregion_phantom_plan.md`;
   open dual-head (`phi_gaps`) and mixed BGO/CsI (per-block crystals) geometries — engine gates.
4. Docstrings + Literate/Documenter doc-site; `load_xcom` hardening + CsI(Tl) config.
