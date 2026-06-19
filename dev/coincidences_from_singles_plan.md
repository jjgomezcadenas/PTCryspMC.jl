# Plan — coincidence selection core + `build_coincidences_from_singles.jl`

**STATUS: DONE.** `src/coincidences.jl` (shared core) + `src/coincidences_hdf5.jl` (LOR HDF5) +
`build_coincidences.jl` (minimal edit, byte-identical) + `scripts/build_coincidences_from_singles.jl`
(singles either-format → `lors_{truth,det}.h5`) + `scripts/tests/check_lors.jl`. `Pkg.test` 770
(incl. the singles-fill ≡ full-stack-fill equivalence + LOR HDF5 round-trip). Honest result:
discrete fields exact, energy to float precision (summation order). Deferred: real times + the
`random` LORs (Step 5); a Julia LOR-HDF5 reader/reducer. See `dev/dev_steps.md`.

Self-contained plan (pick up after a context clear). Goal: form LORs from the **singles
stack** (`simulate_source_mt.jl` output), reusing the existing selection. Structure mirrors the
dev/production split we already have:

| stage | dev (full stack) | production (singles) |
|---|---|---|
| transport | `simulate_phantom.jl` → `stack.csv` | `simulate_source_mt.jl` → `singles.{csv,h5}` |
| LORs | `build_coincidences.jl` → `coincidences_*.csv` | **`build_coincidences_from_singles.jl` → `coincidences_*.h5`** (new) |

Both LOR builders share one selection core in `src/`.

## 1. `src/coincidences.jl` — the shared selection core

Lift the selection out of `build_coincidences.jl` (currently script-local) into the package,
**sink-agnostic** (emits each accepted coincidence via a callback — script → CSV/HDF5 row,
test → vector), so both drivers and the tests run the same code.

- `Response(sigma_xyz, eres, emin, apply_window, win_half)` + `pass_energy(e, r)` — moved as is.
- `GammaAcc` + `reset!` + `contained_one` — moved as is.
- `fill_full!(acc, x, y, z, e_dep, iz, iphi, phscat_row)` — accumulate one full-stack scanner
  deposit (first deposit = LOR point, sum energy, flag overspill). The current inline logic.
- `fill_singles!(acc, x, y, z, e, iz, iphi, nblocks, phscat)` — set `GammaAcc` directly from a
  singles row (`reached=true`, `overspill = nblocks>1`). The singles row *is* the hit.
- `finish_event!(emit, ev, g1, g2, x0, y0, z0, resp, rng) -> (emitted, is_true)` — smear (σ_xyz,
  FWHM(E)) + energy window + two-block selection; on accept calls `emit(ev, hit1…, hit2…, x0,y0,z0,
  truth)`. No I/O. The single source of truth for the selection physics.

Coincidence HDF5 output (parallel to singles), in `src/coincidences.jl`:
- `CoincidenceBuffer` (columnar, quantized): `event` Int32, `truth` Int8 (0=true,1=scatter,
  2=random), per hit `x,y,z` Int16 @0.1 mm + `e` Int16 @0.1 keV + `iz,iphi` Int16 + `t` Float32
  (dummy 0.0; real ns encoding is a Step-5 decision, deferred), emission `x0,y0,z0` Int16.
  Energies are **smeared** here → encode with **clamp** to [0, range] (a Gaussian tail can push a
  small energy negative), unlike the singles truth-energy encode which bounds-*errors*.
- A streaming `CoincidenceWriter` over extensible HDF5 datasets (the total is unknown until the
  stream ends): buffer ~2²⁰ rows → append a block → set_extent. Chunked + shuffle+deflate-4,
  with run params as root attributes. Bounded memory.

## 2. `build_coincidences.jl` — minimal edit (behavior unchanged)

Delete its local `Response`/`GammaAcc`/`pass_energy`/`contained_one`/`finish_event!`; use the
`src` versions. Its full-stack streaming loop now calls `fill_full!` and passes a **CSV-row
sink** to `finish_event!`. Output schema and values **identical**.
**Regression proof:** run it on an existing `stack.csv` before and after, `diff` the
`coincidences_*.csv` → byte-identical.

## 3. `build_coincidences_from_singles.jl` — new driver

Singles → coincidences, single-threaded streaming pass.
- **Input** singles in **either format** (default `output/<tag>/singles.{h5,csv}` by
  `[output].format`; `--singles` override; format by extension).
  - HDF5: `foreach_singles_hdf5` batches → decode Int16→mm/keV → `fill_singles!`; the
    event-grouping state (`cur_ev`, `g[1]`, `g[2]`, `x0…`) persists across batches.
  - CSV: line parser over the singles schema → `fill_singles!`.
- **Selection**: `Response` from `[detector]`; `finish_event!` per event with a sink that
  `push_coincidence!`s into the `CoincidenceWriter`.
- **Output** `output/<tag>/lors_{truth|det}.h5` (truth mode if all detector knobs off). Named
  "lors" because the deliverable is the list-mode LOR list (one record = one LOR). `truth ∈
  {true=0, scatter=1}` now; `random=2` is reserved for Step 5. A root attribute
  **`has_randoms=false`** flags that this is the LOR list *minus randoms* (the complete,
  time-ordered measurement = these LORs ∪ Step-5 random LORs, merged once real times exist).
  The dev `build_coincidences.jl` keeps its `coincidences_{truth,det}.csv` name (plotters).
- **Acceptance denominator**: total `nevents` from the HDF5 `nevents` attribute (or config) —
  the singles stack has no rows for both-missed events, so "events seen" undercounts.

(The dev CSV coincidences are still read by the Python plotters. Reading/reducing the
*production* coincidence HDF5 should be **Julia** — typed, fast, reusing the `src` HDF5
machinery, same reason `build_coincidences` is Julia — emitting small reduced arrays; Python
only for final plotting of that reduced data, if at all.)

## 4. Tests (`test/runtests.jl`)

- **Shared-core equivalence (the payoff):** a handful of events → build a full stack
  (`navigate_photon` → `fill_full!`) and its faithful singles reduction (`_reduce_steps` →
  `fill_singles!`) for the *same* events; run the selection core on each; assert **identical
  emitted coincidences** (truth mode exact; det mode exact too — same `GammaAcc` → same smear
  RNG). The bit-for-bit check we deferred.
- **Coincidence HDF5 round-trip:** push coincidences → write → read back → compare.

## 5. Validation (integration)

- `build_coincidences.jl` byte-diff before/after the minimal edit.
- Run `build_coincidences_from_singles.jl` on a real `singles.h5`; sanity (acceptance %, truth
  split) and a `check`/read of the coincidence HDF5.

## Files

- `src/coincidences.jl` (new) — selection core + `CoincidenceBuffer`/`CoincidenceWriter` (HDF5).
- `src/PTCryspMC.jl` — include + exports.
- `scripts/build_coincidences.jl` — minimal edit to use the core.
- `scripts/build_coincidences_from_singles.jl` (new) — the singles → coincidences driver.
- `test/runtests.jl` — equivalence + HDF5 round-trip testsets.
- `dev/dev_steps.md` — log it.

## Deferred

- Real **times** (Step 5 randoms) — `t` stays Float32 0.0; the `random` truth code (2) is
  reserved for the randoms pass.
- A **Julia** reader/reducer for the production coincidence HDF5 (analysis/reconstruction
  side; faster than Python). Python kept only for final plotting of already-reduced data.
- CSV output option for the singles driver (add only if dev-plotting singles-derived
  coincidences is wanted).
