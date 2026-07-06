# API source scenario — implementation plan

**✅ BUILT & VALIDATED — historical build record.** All 8 steps below are done and merged to `main`
(the `api-scenario` branch is gone). This file is kept as the design rationale (escaped-positron
handling, prompt-gamma deferral, the phase-1/phase-2 seed split); for the current state see
`dev/status.md`. Each step's ✅ note records what shipped.

Builds the **Proton Activity (API)** source mode: read a frozen
`ptcryspg4` scenario (emitters + per-isotope decay budget + phantom) and drive the existing
photon engine from it. Target scenario: `uniform_headep_sobp_1e8` (single-region brain
ellipsoid). Design context in `docs/PTCryspMC_app.tex` §3 and `dev/status.md` "API build plan".

**Decisions (locked):** phantom built from `phantom_regions.csv` (the scenario is the single
source of truth); scanner axis along the beam (+z); Poisson draw done natively in Julia (not a
port of `budget_gen.py` — that stays a statistical reference only, RNG differs); `xcom_brain.csv`
already fetched from NIST into `data/`.

**Scenario path:** `~/Projects/ptcrysp-scenarios/scenarios/uniform_headep_sobp_1e8/`.

---

## Scope boundaries (deferred, with a clear error where hit)

- **Multi-region phantom** (`headep`, `mird_head`): the navigator has one phantom leaf. This
  scenario is single-region, so the loader reads `phantom_regions.csv` generically but **errors
  if >1 region** ("multi-region phantom deferred"). The generic read keeps the single-source-of-
  truth and marks the extension point. **Fully scoped for a future instance in
  `dev/multiregion_phantom_plan.md`** (data, the two extra materials + XCOM-table prerequisite, the
  `MultiRegionPhantom` + region-aware navigator design, tests) — bone attenuates ~1.7× brain, the
  reason it matters. Not on the critical path for the uniform-phantom API source below.
- **Region rotation** (Euler angles): 0 in every current scenario. Loader **errors if any Euler
  angle ≠ 0**. Matches the existing deferred-rotation nit.
- **Prompt-gamma isotopes** (C10, O14 carry `prompt_gamma=1`): the de-excitation γ in coincidence
  with the decay is ignored on the first pass (~2% of the fast budget). Stamped in provenance as
  `prompt_gamma_modeled=false`.
- **Escaped positrons — DROPPED entirely** (decision). ~0.38% of emitters annihilate *outside*
  the phantom: 0.37% surface-pinned (positrons that would leave into air, clamped to the boundary
  by ptcryspg4 — clustered at the beam entrance face) + 0.009% air-produced. A positron ranges
  ~800× further in air (ρ⁻¹), so it annihilates tens of cm away (across the bore / in the crystal /
  out the ends), *not* at the pinned point — these are lost positrons, not a boundary source, and
  keeping them would paint a spurious ~0.4% activity shell on the phantom surface. **Drop criterion:
  `is_inside(phantom, anh)` false** (reuses the `Ellipsoid`). Applied in the scenario reader (step 4)
  by filtering the per-isotope pools at load; the effective count uses the per-isotope inside
  fraction, `M_j ~ Poisson(N_j · f_inside_j)`, so the loss is modeled (not silently over-counting
  survivors). Config flag `keep_escaped` (default false) for a one-off sensitivity check; dropped
  count stamped in provenance.
- **Direction-level common-mode independent of nchunks**: not needed for the first branch (pinned
  nchunks already gives it). A per-event counter RNG is the future strengthening.

---

## Step 1 — `Ellipsoid <: Solid`  (engine, isolated)

`src/geometry.jl`. Mirrors `Sphere`, three semi-axes.

```
struct Ellipsoid <: Solid
    a_cm::Float64; b_cm::Float64; c_cm::Float64   # semi-axes along x,y,z (local frame, origin-centred)
end
volume(e)    = (4/3)π·a·b·c
is_inside(e,p) = (p1/a)^2 + (p2/b)^2 + (p3/c)^2 <= 1
```

Ray crossings (mirror `_sphere_crossings`): substitute `u_i = pos_i/s_i`, `v_i = dir_i/s_i`
into `Σ((pos_i + t·dir_i)/s_i)^2 = 1` → quadratic in the **physical** ray parameter `t`
(dir is unit, t stays cm): `A=Σv_i²`, `B=2Σu_iv_i`, `C=Σu_i²−1`. Same root selection as the
sphere (t > SURFACE_EPS; near/far). `distance_to_exit/entry` identical structure.

`load_solid`: add `shape == "ellipsoid"` → `Ellipsoid(a_cm,b_cm,c_cm)`; extend the error list.

**Tests** (`test/runtests.jl`, new testset mirroring "sphere solid"): is_inside on/off axis;
a ray along each axis exits at the matching semi-axis; an off-axis chord; a miss; entry from
outside; degenerate a=b=c reproduces the `Sphere` results (cross-check). Reduce to `Sphere`
when axes equal.

---

## Step 2 — Brain material  (data, trivial)

`data/materials.json`: add
```
"Brain": { "density": 1.04, "components": {H:.107,C:.145,N:.022,O:.712,Na:.002,P:.004,S:.002,Cl:.003,K:.003},
           "xcom": "xcom_brain.csv" }
```
`xcom_brain.csv` is present (NIST, water-format). No scintillator props (phantom material).

**Validation test:** `sigma_macro(Brain, 0.511)` total (incoherent+photoel+pair) × ... vs the
scenario meta μ = 0.09913 cm⁻¹ (μ/ρ 0.09532 × ρ 1.04). The transport uses incoherent+photoel+pair
(no coherent), so compare against the matching channel sum; expect agreement to the table grid.

---

## Step 3 — Phantom from `phantom_regions.csv`  (loader)

New reader (`src/source.jl` or a small `src/scenario.jl`): parse `phantom_regions.csv`, mm→cm.
Single region → a `PhysicalVolume{Ellipsoid}` (or Cylinder, per `solid` col) at `(cx,cy,cz)/10`
with the region's material. Errors on >1 region and on any nonzero Euler angle (see boundaries).
The world (Air mother) and scanner still come from a geometry JSON; only the **phantom** section
is replaced by the scenario-derived volume. So the API driver composes: world+scanner from a
`geometry_head.json`, phantom from `phantom_regions.csv`.

**Tests:** load the scenario's regions file → an Ellipsoid phantom with the right axes/centre/
material; a synthetic 2-region file throws; a nonzero-Euler file throws.

---

## Step 4 — Scenario reader  (data → source inputs)

`src/scenario.jl`: a `Scenario` struct holding the parsed handoff —
- per-isotope emitter pools: read `emitters.csv` once, group `anh_{x,y,z}_mm` by `isotope_id`
  (mm→cm) → `Vector{Vector{NTuple{3,Float64}}}` (or SoA); **filter to `is_inside(phantom, anh)`**
  unless `keep_escaped` (drop the escaped positrons, see boundaries), recording `f_inside_j` and the
  dropped count per isotope;
- `isotopes.csv` → id → (name, half_life_s) into the `Isotope`/`ISOTOPES` machinery;
- `sampling_budget_<budget>.csv` → `N_expected[id]`; `_meta.csv` → `t_meas_s`, `target_dose_Gy`;
- `run_meta.csv` → provenance (geometry, phantom_material, seed, geant4_version …).

Stamp `scenario` name + `budget` into every output attr. Dose rescale: `N_j(D) = N_expected_j ·
D / target_dose` (linear).

**Tests:** load the real scenario → pool sizes match `wc -l` per isotope; N_expected matches the
CSV; dose rescale linear.

---

## Step 5 — `APISource` + phase-1 materialization  (source)

`src/source.jl`:
```
struct APISource <: Source
    points::Vector{NTuple{3,Float64}}   # the drawn annihilation points, cm  (event_number → point)
    isotope::Vector{Int8}               # per-event isotope id, parallel to points
    lambdas::NTuple{Niso,Float64}       # per-isotope decay constant (for event_time)
end
```
**Phase 1 (build):** given the scenario, budget, `master_seed`, `realization` —
`rng_src = MersenneTwister(master_seed + realization)`; for `j` in id order:
`M_j = rand(rng_src, Poisson(N_j))` (native Julia — Poisson via `Distributions` or an inline
Knuth/transformed sampler to avoid the dep — decide at build; N_j ~ 5e7 so a normal-rounded
draw is acceptable and dependency-free), then `M_j` indices uniform-with-replacement into pool j
→ append point+isotope. `event_number = 1..N` in construction order. Returns the `APISource`.
Depends only on `(emitters, budget, master_seed, realization)` — nchunks/detector-independent.

**Emit split** (so the count-driven/clinic paths are untouched): factor `emit_pair` into
`sample_position` (drawn) + direction draw. Add `event_point(src, ev, rng)` dispatch —
`sample_position(src,rng)` for drawn sources (ignores ev), `src.points[ev]` for `APISource`
(ignores rng) — and `event_isotope(src, ev)` similarly. `singles_chunk!` calls these, then draws
the two back-to-back directions from the chunk rng. One loop, alloc-free, dispatch on source type.

**Tests:** phase-1 with a fixed (master_seed, realization) is reproducible and nchunks-invariant
(build with the array, then verify chunk_ranges over it is independent of nchunks for the POINTS);
`ΣM_j` ≈ ΣN_j within Poisson tolerance; per-isotope M_j mean/var ≈ N_j over many seeds.

---

## Step 6 — Isotope singles column + per-isotope timing  (schema + randoms)

- ✅ `isotope` (Int8) added to the singles schema (`src/singles.jl` `singles_columns`/`singles_doc`,
  binary part I/O, `src/singles_hdf5.jl`), exactly as `n_scatter` was; populated for all modes via
  `event_isotope(src, ev)` (0 for the drawn count/clinic sources). CSV header + `simulate_source_mt`
  sinks updated. `docs/SCHEMA.md` regenerated (drift test green).
- ✅ Per-isotope `scenario_activity_models(scn; seed)` → a `Vector{ActivityModel}` (window [0,t_meas],
  λ from each half-life). Verified: per-isotope truncated-exponential mean matches the analytic to
  <1% (C10 front-loads to ~28 s, C11 to ~530 s in a 1200 s window). `scripts/tests/check_scenario.jl`
  reports + gates it.
- **DEFERRED to step 7** (needs the driver's scenario models): `build_randoms` reading the isotope
  column to time each single by its isotope's λ. The column and the model helper are in place; the
  wiring (API-mode model construction + the randoms per-isotope restore) lands with the driver.

**Tests (done):** singles round-trip carries isotope; `singles_chunk!` stamps each single's isotope
via `event_isotope`; `scenario_activity_models` gives the right per-isotope truncated exponential.

---

## Step 7 — `mode="api"` driver branch + config  (driver)  ✅ DONE

Built and verified end-to-end. `simulate_source_mt` has a `mode="api"` branch: load scenario →
phantom from the scenario (overrides the JSON placeholder) → materialize the `APISource`
(master_seed+realization, separate from the transport seed) → transport. The singles carry the
`isotope` column and the run stamps `source_mode/scenario/budget/dose/realization/master_seed` +
the per-isotope timing (`isotope_half_lives`, `isotope_names`, `t_meas_s`, `time_seed`) into the
HDF5 attrs. `build_randoms` reads those attrs → per-isotope `ActivityModel`s, times each single by
its isotope's λ (falls back to the single config model when the attr is absent — clinic/count
unchanged). `runs/uniform_headep_bgo_api.toml` + `geometry/geometry_head.json` (world+ring; the
phantom section is a placeholder the driver replaces). `run_prod.sh` skips `--nevents` for `api`
(N is materialized). Verified: the full chain (simulate → truth → randoms → reco) runs and
`check_lors` PASSes on the API `lors_det.h5`; the isotope column populates in the budget
proportions (O15:C11 ≈ 2.86 vs 2.89 expected); randoms 2τS² validation at full dose is step 8.

---

## Step 8 — Validation  ✅ DONE

Reference run: `uniform_headep_sobp_1e8`, budget `fast`, dose 1 Gy, BGO, master_seed 1 / realization 0
→ 80.18M events, the full chain in ~5 min on 18 threads. `scripts/tests/check_api_validation.jl`
gates it, all green:
- **N & per-isotope**: materialized N = 80,177,205 = singles `nevents`; every M_j within ~2σ of its
  Poisson budget.
- **Reproducibility / nchunks-invariance**: every singles emission point equals the materialized
  array's `point[event_number]` to the 0.1 mm quantization — the source is indexed by event, not
  redrawn per chunk.
- **Spatial activity profile**: source↔detected z-profile correlation = 1.0000 (detection does not
  distort the SOBP activity shape); distal 50% activity edge z ≈ −18 mm (proximal to the target
  distal face z = −5 mm, the physical activity-vs-dose offset the range study measures).
- **Randoms 2τS² with per-isotope timing**: 84,820 measured vs 85,173 analytic → **ratio 1.00**
  (the new multi-isotope timing reproduces the accidental-coincidence law at full dose).
- **LOR split**: 46.37M LORs — 32.1% true / 67.7% scatter / 0.18% random (physical for a brain head).

The reference `prod/uniform_headep_bgo_api/lors_det.h5` is the first end-to-end API deliverable —
ready for the downstream reconstruction / R study.

---

## Build order & dependencies

1 (Ellipsoid) and 2 (Brain) are independent, testable in isolation → do first, together.
3 (phantom loader) needs 1. 4 (scenario reader) is independent data. 5 (APISource) needs 4.
6 (isotope column) is schema, needs the source to carry isotope (5). 7 (driver) composes 3+5+6.
8 validates end-to-end. Each step lands with its tests green before the next.
