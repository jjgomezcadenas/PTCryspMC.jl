# Project status — PTCryspMC.jl

A one-page snapshot of where the simulation stands, plus the **deferred-work register**.
Companion docs: the *method* is in `docs/PTCryspMC_phys.tex` (engine) and `docs/PTCryspMC_app.tex` (modes); the
decisions + code layout in `CLAUDE.md`.

_Last updated: 2026-06-26._

## What it does

A fast, photon-only Monte Carlo: from a **source** (a clinical tracer distribution, or — planned —
the positron activity a proton field leaves) + a detector description, write the list-mode
coincidence/LOR file a PET scanner would record. It never runs proton transport.

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
  (`jitter = −ln u / (N_det·r0)`, `N_det = yield·E·pde`, `r0 = Σ wₖ/τₖ`), stamped **once** at
  generation, **relative to the decay** so it stays small (Float32). Absolute time =
  `event_time(ev)·1e9 + t_rel`, reconstructed only where randoms need a common clock.
- **Coincidence window τ = 3 ns** (CsI and BGO), read off the DT study (`examine_dt.jl` +
  `py/plot_dt.py`); `compare_crystal_timing.jl` explains why τ is crystal-independent.
- **Reco lower energy cut** `reco_emin_keV = 450` (photopeak region; the spectrum studies keep
  `emin_keV = 300` to see the Compton shoulder). No upper cut (it's an analysis-time choice).
- **Validated:** `Pkg.test` **870**; randoms match the analytic `2τS²` (CsI 1248 vs 1291, ratio 0.97;
  the clinical Vacuum/BGO 10⁸ run 52694 vs 52708, ratio 1.00); reco acceptance CsI 8.98% / BGO 23.75%;
  corrected residual median|dt| (CsI 0.059 / BGO 0.200 ns) matches the analytic single-photon jitter;
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
  `sphere_water_f18_csi`, `nema_iq_f18_csi`, `sphere_air_bgo`. `run_prod.sh` derives N for clinic
  configs (no `--nevents`).
- **API** (After Proton Irradiation, count-driven) — *planned* (the second source branch). A frozen
  `ptcryspg4` scenario supplies the emitters and the given per-isotope decay budget. Design in
  `docs/PTCryspMC_app.tex`; build plan below.

### API build plan (the next feature — continuity notes)

Parallels Clinical (same source abstraction → N annihilations, per-event point, per-event isotope→λ),
but count-driven and multi-isotope. **Blocked on the upstream `ptcrysp-scenarios` format** — the
reader columns can't be fixed until it is (see `CLAUDE.md` "Input — a scenario"). Build order:

1. **Scenario reader** (format-pending): `emitters.csv` (annihilation point + isotope tag, positron
   range already applied), `sampling_budget_<scn>.csv` (`N_j` per isotope, given), `isotopes.csv`
   (code → `T½`, `β⁺`), `run_meta.csv` (normalization). Stamp the scenario name into every output.
2. **Per-scenario isotopes**: feed `isotopes.csv` into the `Isotope` table (already `(T½, β⁺)`, in
   `src/activity.jl`).
3. **`APISource <: Source`** (in `src/source.jl`, mirroring `ClinicSource`): holds the emitters +
   budgets; `sample_position` draws an emitter (with replacement) per isotope to `N_j`, returning the
   point and its isotope; `N = Σ N_j` (apply `β⁺` per isotope if the `N_j` are nuclear decays —
   format-dependent).
4. **Per-isotope timing + an isotope column on the singles** (a singles-schema add, exactly like
   `n_scatter` was): `simulate_source_mt` records each single's isotope; `event_time`/`ActivityModel`
   generalize to use the emitter's isotope λ; `build_randoms` pairs cross-decay singles with the
   per-single λ. `SCHEMA.md` regenerates from the code.
5. **Driver**: a `[source].mode="api"` branch in `simulate_source_mt` (load scenario → `APISource` →
   per-isotope `ActivityModel`), alongside the existing clinic/count-driven branches.
6. **Validate** against a real scenario; reuse the seed/`nchunks` reproducibility (the emitter draw
   must be seeded per chunk).

Buildable now regardless of format: the per-isotope `ActivityModel`, the singles `isotope` column, the
`APISource` shape. Format-blocked: the reader (1–2) and whether `N_j` carry `β⁺`.

## Detector configs

CsI (5% FWHM @ 511 keV) and cryogenic BGO (10%) — sphere-water scenarios, validated. CsI(Tl) (7%)
and the pixelated detectors (LYSO, standard BGO) are still to add.

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

- **API source scenario** — the second source branch (Clinical is built; see "Source scenarios").
  `[source].mode="api"`, `scenario=…` → a frozen `ptcryspg4` scenario: `emitters.csv` (points +
  isotope, range pre-applied), the measured per-isotope budget `N_j` (given, not derived),
  `isotopes.csv` (per-isotope `T½`). Needs the scenario reader; adds an isotope column to the singles
  + per-isotope `event_time`. Multi-region clinic spatial draw uses dynamic dispatch per event — worth
  a glance if a clinic run ever goes to 10⁸ with many regions.
- Refine the per-crystal **PDE** (0.45 placeholder for CsI and BGO; should differ by emission colour).
- Cylinder/vacuum configs need their own `examine_dt` (different scanner length → different TOF tail) before fixing their `[timing].tau_ns`.
- Threshold / CFD timing (the "first photon" model is the leading-edge idealization).
- Pixelated detectors report a fixed crystal, not a continuous position. Placement **rotation** transform (only when a volume needs it).

## Next

1. **API source scenario** — the second source branch (the `ptcryspg4` scenario reader), when the
   upstream format is fixed.
2. **Docstrings + Literate/Documenter doc-site** (the remaining documentation task).
3. Crystal XCOM tables + the `load_xcom` hardening; the CsI(Tl) config.
4. Downstream, separate & deferred: range precision, detector comparison, reconstruction.
