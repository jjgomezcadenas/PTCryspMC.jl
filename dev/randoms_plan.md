# Plan — Step 5: randoms (the third LOR category) + the timing model

> **Status: implemented.** Step 5 is done (timing, randoms, reco) — see `dev/status.md` for the
> current state. Divergences from this plan: the truth builder is `build_true_coincidences_from_singles.jl`
> (truth-only); `t_rel` is stamped once in the singles, not in the LOR builder; `dt` is the signed
> `(t1−t2) − TOF_diff`. This file records the original design.

Self-contained (pick up after a context clear). State going in: the production chain
`simulate_source_mt.jl` (MT) → `prod/<tag>/singles.h5` → `build_coincidences_from_singles.jl`
→ `prod/<tag>/lors_{truth,det}.h5` is done (trues + scatters; `truth ∈ {0,1}`, `random=2`
reserved, `has_randoms=false`, `t1/t2` dummy `0.0`). **Piece 1 (the activity model,
`src/activity.jl`) is done.** Randoms are the only missing LOR category and complete the
list-mode measurement.

## What a random is

Two singles from **different** annihilations whose timestamps fall within the coincidence
window τ → a false LOR. Rate ≈ **2·τ·S²** (S = singles rate), **front-loaded** (∝ activity²),
a few % of trues. Method: time-stamp each single, sort by time, pair **cross-annihilation**
singles within τ → `truth=2`.

## Timing model — the first detected photon (NOT a Gaussian σ_t)

A gamma's timestamp is the arrival of its **first detected scintillation photon**, a physical
distribution — not a free σ_t. For a deposit of energy `E`:

- detected photons `N_det = light_yield · E · PDE` (e.g. cryo CsI ~10⁵ γ/MeV × ~0.5 MeV × PDE
  → ~10⁴);
- each scintillation photon is emitted ~`Exp(τ_scint)` after the deposit (τ_scint ≈ 1 µs);
- the **first** detected photon = min of `N_det` such draws → `Exp(N_det/τ_scint)`, i.e.
  `jitter = −(τ_scint/N_det)·ln u`, mean `τ_scint/N_det` (~0.1 ns here). The per-gamma timing
  jitter falls out of yield/PDE/τ_scint/E — no tunable σ_t.

So per single (gamma `g`):
```
t0(g) = t_annih(ev) + TOF(g) + jitter(E_g)
        t_annih = event_time(model, ev)              # the activity model (piece 1)
        TOF(g)  = ‖hit_g − emit‖ / c                  # c ≈ 299.792 mm/ns
        jitter  = −(τ_scint / N_det) · ln u           # N_det = yield · E_g · PDE
```

### DT, TOF, and the window (the user's scheme)

For a **true** pair (same annihilation, common emit `(x0,y0,z0)`):
```
TOF_diff = ( ‖hit1 − emit‖ − ‖hit2 − emit‖ ) / c     # geometric, from fields already stored
DT       = |t0g1 − t0g2| − TOF_diff                  # → 0 with no fluctuations; this IS the
                                                     #   timing resolution (the jitter residual)
```
**TOF is required to form DT** (subtracts the geometric spread), and it's computable per pair
from `(x0,y0,z0)` + the two hits the LOR already carries. The **coincidence window τ** is then
set from the **trues' `|Δt0|` distribution** — chosen for high true-pair acceptance — and that
same window, applied to the singles stream, lets the accidental cross-event singles **leak in
as randoms**. (We don't separately reason about FOV-TOF vs resolution; we read τ off the trues.
Randoms have no common emit → no TOF_diff/DT → they're paired on raw `|Δt0|` within τ.)

## Locked decisions

- **`event_time(m, ev)`** — per-event annihilation time, deterministic, alloc-free, keyed by
  `event_number` → MT-chunking-independent (piece 1, done).
- **Timestamp = scintillation first-photon model** (above), in a new `src/timing.jl` — replaces
  the Gaussian `smear_time` idea. Params: `light_yield` [γ/MeV], `tau_scint_ns`, `PDE`.
- **DT needs TOF**, computed per pair from `(x0, hit1, hit2)`. Window τ derived from the trues.
- **Two scripts:** the trues/scatters builder stamps `t1,t2,dt`; a new single-threaded
  `build_randoms_from_singles.jl` does the cross-annihilation pairing (global time-sort).
- **Build single-threaded first**, benchmark on the 10⁷ file, MT only what's justified.
- **Complete measurement = trues+scatters ∪ randoms, time-ordered**, one `lors.h5`,
  `has_randoms=true`.

## LOR schema additions (per pair, one row)

Populate/add to the existing `CoincidenceBuffer` (`src/coincidences_hdf5.jl`):
- **`t1_ns`, `t2_ns`** — already in the schema (dummy 0.0); now the two timestamps `t0g1,t0g2`.
- **`dt_ns`** — new per-pair column, `DT` (the resolution residual). `NaN` for randoms.
- **`tof_diff_ns`** — optional per-pair column; fully derivable from `(x0,hits)` so skippable,
  but keep it if we want the file self-contained. `NaN` for randoms.
Pair quantities are stored **once** per row (no `_1`/`_2` repetition). Unit suffixes (`_ns`)
match the schema convention (`x_mm`, `e_keV`).

## Pieces (ordered; ST first)

1. **`src/activity.jl`** — DONE (toy ¹⁵O activity + `event_time`).
2. **`src/timing.jl`** — pure functions; the **crystal `Material` already carries** `light_yield`,
   `scint_decay_ns`, `scint_decay_w` (DONE), so the only extra input is **PDE** (readout, config).
   - `first_photon_jitter(mat, E_MeV, rng)` [ns]: `N_det = mat.light_yield·E_MeV·mat.pde` (PDE is
     per-crystal in the DB too); the first of `N_det` photons from the decay mixture. The first
     photon lands at ~0.1 ns ≪ τ, so it's the min of `N_det` exponentials at the **effective
     initial rate** `r0 = Σ wₖ/τₖ` → `jitter = −ln(u)/(N_det·r0)` (exact in the early-time limit;
     handles 1 *or* 2 components uniformly). Mean `1/(N_det·r0)`.
   - `tof_ns(emit, hit)` = `‖hit−emit‖/c`, `c = 299.792458 mm/ns`.
   - `photon_timestamp(t_annih, emit, hit, E_MeV, mat, rng)` = sum of the three.
   Tested: mean jitter ≈ `1/(N_det·r0)` for CsI (single) + BGO (two-component); alloc-free.
3. **Trues/scatters: stamp `t1,t2,dt` + report the window.** Extend
   `build_coincidences_from_singles.jl`: time-stamp each gamma, compute `TOF_diff`/`DT`, write
   `t1_ns,t2_ns,dt_ns` (+ optional `tof_diff_ns`). Print the trues' `|Δt0|`/`DT` distribution to
   guide τ (or auto-suggest τ at a high-acceptance quantile).
4. **`scripts/build_randoms_from_singles.jl`** (single-threaded) — keep good singles
   (contained-one + energy selection), time-stamp each, **sort by time**, scan the τ-window
   pairing **cross-event** singles → random LORs (`truth=2`, `dt_ns=NaN`) → `randoms_{mode}.h5`.
   Multiples policy (reject >2-in-window vs naive pairwise) decided at impl.
5. **Merge → complete `lors.h5`** — trues+scatters ∪ randoms, **sort by coincidence time**,
   `has_randoms=true`. A `merge_lors.jl` or appended by `build_randoms`.

## The coincidence window τ

Set from the **trues' `|Δt0|`** (piece 3 reports it) for high true-pair acceptance; a config
`[timing].tau_ns` carries the chosen value (or piece 3 auto-suggests, e.g. the quantile keeping
99 % of trues). The same τ governs the randoms pairing.

## Pairing at scale

10⁷ → ~13.7 M singles. First cut: **in-memory global sort** of the good singles' times. For
10⁸: **τ-bin streaming** (bin by time, pair within/adjacent bins) — the natural MT partition.
Cross-event only; no opposition filter.

## Benchmark (10⁷ file) → MT decision

- `build_coincidences_from_singles` (~5 s ST): split read/decompress vs select+timestamp vs
  write; MT (event-range chunks, reuse `simulate_source_mt` pattern) only if CPU-bound.
- `build_randoms` sort+pair: prototype ST, extrapolate to 10⁸; sort is the likeliest MT win
  (τ-bins / parallel sort).
- Deciding ratio: each step's ST 10⁸ time vs the transport's (~30–60 s); under ~10 s → skip MT.

## Validation

- **Random rate** ≈ analytic `2τS²` (vary τ → linear; vary S via `t_acq` → quadratic).
- **Front-loaded** — randoms-vs-time ∝ activity² (faster than trues ∝ activity).
- **DT distribution** = the timing resolution; its width ≈ `√2 / (N_det·r0)` (two gammas,
  `r0 = Σ wₖ/τₖ`), scaling as expected with yield/PDE/E. An example study + plot (à la `o15_lifetime`).
- Extend `scripts/tests/check_lors.jl` to the three-way split (true/scatter/**random**), the DT
  spread, and `has_randoms=true`.

## Config additions

```toml
[timing]   t0_s, t1_s, half_life_s, time_seed, tau_ns   # activity window + coincidence window
```
No `[detector]` timing cards: `light_yield_per_MeV`, `scint_decay_ns`/`scint_decay_w`, `eres_a`,
**and `pde`** are all keyed to the crystal in `data/materials.json` (DONE). PDE is per-crystal
because it's the photodetector efficiency *at the crystal's emission wavelength* — currently a
**0.45 placeholder for both CsI and BGO; should differ by emission colour (CsI UV vs BGO ~480 nm)
— refine with real numbers.**

## Deferred

- **Real per-isotope activity** — Step 1 (scenario source: isotope tags + budget). Swapping it
  in changes only `event_time` + adds an isotope column to the singles.
- **Threshold / N-th-photon / CFD** timing — "first photon" is the leading-edge idealization;
  a discriminator model is a small extension.
- **MT** for whichever consumer the benchmark justifies.
- A Julia reader/reducer for the production LOR HDF5 (analysis side).
