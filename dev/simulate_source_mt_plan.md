# Plan — `simulate_source_mt.jl` (production singles pipeline)

**STATUS: DONE** (2026-06-19). `navigate_single_photons` (alloc-free, `src/navigator.jl`) +
`scripts/simulate_source_mt.jl` (MT, singles-only) built and validated: `Pkg.test` 710
assertions pass (incl. zero-alloc + reduction-match); `-t 1` ≡ `-t 18` byte-identical at fixed
`nchunks`; 6.1× on 18 cores at 3 M events. Deferred (separate tasks): the `build_coincidences`
singles-reader (validation #3, the bit-for-bit diff) and HDF5. See `dev/dev_steps.md`.

Self-contained plan so this can be picked up after a context clear.

## Goal

A production driver that scales to ~10⁸ decays: **singles-only** output,
**multi-threaded**, config-driven. `simulate_phantom.jl` stays untouched (the full
per-interaction stack / dev / CNN path). Source stays `UniformVolumeSource`/`PointSource`
— nothing new for spheres/cylinders (the proton-beam scenario only swaps the point sampler
later; same medium, same transport).

## Phase 0 — alloc-free per-photon function (`src/navigator.jl`)

The reason the script can scale. The existing `navigate_photon` builds a `Vector{NavStep}`
per photon (and each `propagate_photon` builds its own `recs` vector) — at 2×10⁸ photons
that is the GC ceiling: stop-the-world GC serializes the threads and `-t 18` stops scaling.

Add a sibling that folds each photon to one summary with **zero heap allocation**, sharing
the same physics kernel (`sigma_macro`, `sample_distance`, `distance_to_exit`,
`sample_interaction` — all pure, already alloc-free). Only the bookkeeping shell differs
(fold into locals instead of `push!` into a vector).

```julia
# navigate_photon ("navigate_photon_full_stack"): records EVERY interaction →
# Vector{NavStep}. Allocates. Used by simulate_phantom.jl.
function navigate_photon(...)

# navigate_single_photons: folds each photon to ONE singles summary —
# (reached, first crystal hit x,y,z + block iz,iphi, summed E, nblocks, phantom_scatter).
# No vector, no closure → zero heap allocation per photon. Used by simulate_source_mt.jl.
function navigate_single_photons(...) -> (reached, x,y,z, iz,iphi, e, nblocks, phscat)
```

- "First crystal hit" = first `:scanner` interaction with `e_dep>0` (the LOR point) — the
  same rule `build_coincidences` uses, so the reduction is faithful.
- `e` = summed energy across scanner deposits; `nblocks` = distinct blocks (1 = contained,
  >1 = overspill); `phscat` = any phantom interaction.
- We did **not** use a callback (one loop, caller decides full-vs-reduced): a closure
  capturing/reassigning the accumulators can box them → per-photon allocation, the exact
  cost we are killing. Two standalone functions sharing the kernel avoid that footgun; the
  bit-for-bit validation (below) guards against the two shells drifting.
- **Verify zero-alloc** with `@allocated` on a small run *before* any MT benchmarking (so we
  never measure GC by accident).

## Singles output — `output/<tag>/singles.csv`

One row per *detected* photon (a miss writes nothing; an overspill photon is still a real
single and is written):

```
event_number, gamma, x_mm, y_mm, z_mm, e_keV, iz, iphi, nblocks, phantom_scatter, x0_mm, y0_mm, z0_mm
```

`x,y,z` / `iz,iphi` = first crystal interaction (LOR point + block); `e_keV` = summed block
energy, **truth/unsmeared** (smearing stays in `build_coincidences`); `nblocks` = distinct
blocks touched (1 = contained, >1 = overspill).

## MT structure

- `nthreads = Threads.nthreads()` (from `-t`, default whatever you launch with).
- `nchunks = 8 * nthreads` (default), **overridable** via `--nchunks`.
- Split events `1:N` into `nchunks` **contiguous** ranges.
- **Pre-allocate RNGs before the loop**: `rngs = [MersenneTwister(seed + (c-1)) for c in
  1:nchunks]` (chunk 1 = base `seed`). Inside the loop threads only *draw* — never construct
  or re-seed an RNG in the hot path.
- `@threads :dynamic for c in 1:nchunks` — fast cores grab more chunks, absorbing the
  Super/Performance core heterogeneity (M5, 18 cores). Each chunk uses `rngs[c]`, streams its
  rows to its own part-file `singles.partC.csv` (no header), keeps a local counter.
- After the loop: write the header to `singles.csv`, concatenate the part-files **in chunk
  order** (→ globally event-ordered, so the existing streaming `build_coincidences` reader
  needs no change), delete the parts. Bounded memory throughout.
- Reproducible over `(seed, nchunks, N)` independent of thread count / scheduling (each chunk
  result depends only on its index-seeded RNG and its event range).

## Config / CLI

Same TOML. Reads `[geometry]`, `[source]`, `[transport]` (incl. `seed`), `[output].dir/tag`.
Ignores `[detector]` (smearing stays in `build_coincidences`). CLI overrides: `--config`,
`--nevents`, `--seed`, `--nchunks`.

## Validation

1. **Zero-alloc:** `@allocated` per photon ≈ 0 for `navigate_single_photons`.
2. **Reproducibility (self-contained):** same `(seed, nchunks)`, `-t 18` vs `-t 1` →
   identical `singles.csv` (diff).
3. **Bit-for-bit vs full stack:** at `--nchunks 1`, chunk-1 RNG = `MersenneTwister(seed)`
   reproduces `simulate_phantom`'s exact draw order → the singles must reduce-match the full
   stack. *Requires `build_coincidences` to read singles — separate task.*
4. **Physics:** aggregate metrics (reached %, spectrum, truth/scatter split) agree with
   `simulate_phantom` within Monte-Carlo error.
5. **Benchmark:** `-t 18` vs `-t 6`. If equal → GC-bound; Phase 0 should have removed that.

## Out of scope (separate / deferred)

- **`build_coincidences` singles-reader** — needed to consume the singles stack + run
  validation #3; its own task.
- **HDF5** — CSV only for now.
- **Isotope tag + sampled time** — Step 5 (randoms), layered on the singles row later.

## Files touched

- `src/navigator.jl` — add `navigate_single_photons` + comments on both; export.
- `src/PTCryspMC.jl` — export.
- `scripts/simulate_source_mt.jl` — new driver.
- `test/runtests.jl` — testset for `navigate_single_photons` (reduction matches a full-stack
  reduce on the same events; zero-alloc check).
- `dev/dev_steps.md` — log the production path.
