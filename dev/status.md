# Project status — PTCryspMC.jl

A one-page snapshot of where the simulation stands, plus the **deferred-work register**.
Companion docs: the *method* is in `docs/PTCryspMC_phys.tex` (engine) and `docs/PTCryspMC_app.tex` (modes); the
decisions + code layout in `CLAUDE.md`.

_Last updated: 2026-07-31._

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
  measured 1.84 ns CTR of Soleti et al. 2024; calibration + references in `latex/scanner_prods.tex`),
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
and raise cost; 3X₀ geometries kept for a thick-crystal study). Only BGO_195K + CsI are
produced (the BGO temperatures are timing-degenerate → 195K is the cheaper representative;
`bgo77k`/`3x0` configs ready but unscheduled). The old hybrid "BGO" material is removed;
`runs/uniform_headep_bgo_api.toml` is the frozen reference of the existing `crysp_ring_1m/bgo`
master. CsI(Tl) (7%) and the pixelated detectors (LYSO, standard BGO) are still to add.

**Scanner-geometry survey — HISTORICAL (off-centre; SUPERSEDED by the v2 generation below; these
masters are DELETED).** A family of closed rings, same block/constants/source, varying radius and
AFOV. Table kept for the geometry-vs-acceptance intuition only (shard 0, off-centre, BGO_195K / CsI):

| scanner | R (mm) | AFOV (mm) | φ×z | accept. BGO / CsI | shards |
|---|---|---|---|---|---|
| ring_1m  | 387 | 1024 | 48×20 | 19.0 / 7.6% | 10 |
| r30_50cm | 300 | 512  | 37×10 | 13.6 / 5.35% | 1 |
| chs      | 200 | 358  | 25×7  | 12.7 / 4.9% | 10 |
| r35_50cm | 350 | 512  | 43×10 | 12.1 / 4.78% | 10 |
| ring_50cm| 387 | 512  | 48×10 | 11.1 / 4.4% | 1 |
| r35_35cm | 350 | 358  | 43×7  | 7.8 / 3.11% | 10 |
| r30_30cm | 300 | 307  | 37×6  | 7.2 / 2.85% | 10 |
| r35_30cm | 350 | 307  | 43×6  | 6.2 / 2.48% | 1 |

Ten full masters (5 geometries × 2 crystals, ΣM = 8.02e8 each) + 3 geometry pairs at shard 0.
Findings: purity geometry-independent (~72–75% / ~87% trues BGO / CsI); radius beats length —
CHS (R200) beats every R≥300 ring below 1 m; the axial-loss hits the phantom poles, not the
central distal edge. `chs` = compact head scanner (R200 → ~30 cm bore, vertex-to-C7).

## Generation-2 (current, 2026-07-13) — supersedes the off-centre survey above

The engine now runs the **v2** products convention (full spec + downstream contract:
`dev/generation2_plan.md`; layout: `dev/PRODUCTS.md`; schema: `docs/SCHEMA.md`):

- **Tumour-centred** patient placement (`[source].center_on="tumour"`, `src/scenario.jl`): the SOBP
  dose-target centre (distal dose R80 − half the target thickness) sits at the ring centre — a fixed,
  dose-based reference, window/mix/washout-independent (+25.6 mm shift here). `"distal_edge"` remains
  in the code but is NOT used (it drifts with the acquisition window).
- **Irradiation-end clock + fixed acquisition scenarios** `[t_del, t_del+t_ac]`, with
  `t_del ∈ {120,180,300}` s and `t_ac = 300` s (config `[timing].t_del_list`/`t_ac_s`). Transport runs
  ONCE over the union window [120,600] s; `reco_lors` band-cuts each scenario → one
  `lors_det_del<NNN>.h5` per scenario. Counts from the irradiation-end population N_j⁰ (budget-
  independent, agrees <0.2%). `t_decay_s` zero moved acquisition-start → **irradiation end**.
- **Isotope truth** per LOR (new `isotope` column) + **self-describing shard metadata**
  (`generation="v2"` guard, embedded `geometry_json`, Mizuno washout `g_i` stamped — NOT applied,
  left to downstream; the full §5 attribute set is drift-tested in SCHEMA). `t_ac=1200` was the
  inherited budget window — replaced by the realistic 300 s.
- `publish_prod`/`run_shards` handle the per-scenario leaves (`del<t_del>s_ac<t_ac>s_<dose>`) + the v2
  pruning. `Pkg.test` green.

**Produced (v2, in PtCryspProds; the tree is now ALL-v2, ~15 GB):** six masters, three scenarios × 10
shards each. **CsI** (r unchanged): `crysp_ring_1m` (R38.7), `crysp_r35_50cm`, `crysp_r35_35cm`.
**BGO_195K** at **R = CsI + 5 cm** (cryostat; new actual-radius geometries): `crysp_ring_1m` (R43.7),
`crysp_r40_50cm`, `crysp_r40_35cm`. Findings (shard 0): acceptance grows with AFOV, falls with delay;
purity geometry/delay-set, ~87% CsI / ~73–77% BGO; **BGO ~2.4× the CsI acceptance** (denser, 413 vs
472 keV cut, net of the +5 cm cryostat). Documented in `latex/scanner_prods.tex`
(fig:csisurvey + fig:crystalcompare, from `py/fig_csi_survey.py` + `py/fig_crystal_compare.py`).
**All old off-centre products (CsI + BGO, incl. the frozen `crysp_ring_1m/bgo` reference) are DELETED.**

**Remaining:** more geometries / the second BGO temperature on demand (`run_shards <v2 cfg> 0 9`);
downstream σ_R (CryspBrainSim vendors `dev/generation2_plan.md` + `docs/SCHEMA.md` + `dev/PRODUCTS.md`
+ `dev/data_generation_strategy.md`, all v2-updated).

## Position resolution model 2 (branch `position-resolution`, 2026-07-17)

The CryspLight two-Gaussian position resolution (`data/pet_resolution_recipe.json`; per-coordinate
core/tail mixture, **averaged over x/y/z** into one isotropic smear) is implemented as **position
model 2** beside the single-Gaussian **model 1** (σ_xyz 1.486 mm), selected per run by
`[detector].pos_model` with an optional selection tier `[detector].selection` = `"none"/"80"/"60"`
(per-gamma Bernoulli efficiency 1.0/0.8/0.6, each tier with its own core fraction; only the
smearing part of the recipe is imported — the efficiency ladder duplicates the transport). Crystal
constants (`pos2_sigma{1,2}_mm`, `pos2_f_core`): CsI-family 1.152/4.226 mm, f = 0.616/0.789/0.865;
BGO (both temps) 0.857/2.560 mm, f = 0.758/0.864/0.891. New **CSI_TL** crystal = CsI(Tl)
(50k ph/MeV, 680/3340 ns 64/36%, 7% eres; CsI attenuation, σ_t 0.35 ns kept, τ 1.5 ns).
**Model 1 remains the default; legacy outputs are bit-identical** (mixture and Bernoulli draw
nothing unless enabled). `Pkg.test` **1039**.

**Produced (model-2 pair, small scanners, 10 shards × 3 scenarios each, ΣM = 4.87e8):**
`crysp_r35_35cm_csi_2x0/csi_tl` — CSI_TL, **80% selection**, window 511±2σ_E (±30.4 keV);
`crysp_r40_35cm_bgo_2x0/bgo_77k` — BGO_77K, **no selection**, window ±43.4 keV (both:
`reco_emin_keV = 0`, the symmetric window is the whole energy selection — first production use of
`window_fwhm`; = 0.84933 = 2σ in FWHM units). Shard-0 numbers (del120): CsI(Tl) acc 1.72%/decay,
purity 89.7%, trues-DCA q50/q68/q95 = 1.74/2.34/4.58 mm; BGO_77K acc 5.35%, purity 86.6%
(vs 76.9% for the lower-cut-only BGO_195K master — the symmetric window), DCA 1.44/1.88/3.29 mm.
Cross-shard acceptance spread ≤0.1% rel. Configs `runs/uniform_headep_r35_35cm_csitl_m2s80.toml`,
`runs/uniform_headep_r40_35cm_bgo77k_m2.toml`. Findings (scratch study + shard 0): model 2 ≈
model 1 for BGO; for CsI the 4.2 mm tail is a real degradation (q95 3.1→5.5 mm no-selection) that
the 80% tier largely buys back at ×0.64 pairs; BGO beats CsI(Tl) on both acceptance and DCA at
matched geometry. **Planned before merging: flip the default to model 2** (JJ 2026-07-17), then
update CLAUDE.md/notes accordingly.

## dd reruns — data-driven source (branch `dd-reruns`, 2026-07-26)

σ_R reruns on the scenario `uniform_headep_sobp_1e8_dd`, whose emitters are sampled from the
**nominal fitted cross sections** (ptcryspg4 phase 3d; record: ptcryspg4
`workshop/xsections_phases.md`). Two configs — `uniform_headep_ring1m_bgo_dd` (TBP reference) and
`uniform_headep_r40_35cm_bgo_dd` (CAFOV) — are the v2 files with two keys changed —
`scenario_dir → …_1e8_dd`, `budget → "d120s300"` (a CBS budget; any budget backs out the same
N_j⁰) — so σ_R and fitted-edge differences from the native masters measure the source alone.
Tumour centring is dose-based and gives the same z offset as native (+25.59 mm): the productions
are registered in z by construction. Scenario pools 1.27M/772k/142k (O15/C11/N13 — the ×1.32 ¹¹C,
×1.49 ¹³N data-driven enrichment); publish under `PtCryspProds/uniform_headep_sobp_1e8_dd/`.

**Timing (this machine, 18 threads, nchunks 144):** one `run_shards.sh` shard of the reference
BGO ring (`uniform_headep_ring1m_bgo_dd`, M = 5.08e7 union-window decays) takes **1 min 53 s**
for the full chain (simulate → coincidences → randoms → reco → publish → gates → prune);
a 10-shard master ≈ **20 min**.

**Produced:** two full masters under `PtCryspProds/uniform_headep_sobp_1e8_dd/`, each 10 shards
(realizations 0–9) × three delay leaves (`del{120,180,300}s_ac300s_1Gy`), gates green, randoms
ratios 0.99–1.01:

- `crysp_ring_1m_bgo_2x0/bgo_195k` — TBP reference (R43.7), ΣM = 5.08e8 per leaf; 9 shards in
  22 min 48 s.
- `crysp_r40_35cm_bgo_2x0/bgo_195k` — CAFOV counterpart (compact BGO, ~35 cm AFOV).

Both feed the downstream data-driven σ_R study (CryspBrainSim `md/results.md`: TBP + CAFOV,
five acquisition protocols).

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
