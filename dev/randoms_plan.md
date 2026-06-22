# Plan — Step 5: randoms (the third LOR category)

Self-contained (pick up after a context clear). Supersedes the pre-singles version of this
file. State going in: the production chain `simulate_source_mt.jl` (MT) → `prod/<tag>/singles.h5`
→ `build_coincidences_from_singles.jl` (single-threaded) → `prod/<tag>/lors_{truth,det}.h5` is
done (trues + scatters; `truth ∈ {0,1}`, `random=2` reserved, `has_randoms=false`). Randoms are
the **only missing LOR category** and complete the list-mode measurement.

## What a random is

Two singles from **different** annihilations arriving within the coincidence window τ (~ns) → a
false LOR with a wrong line of response. Rate ≈ **2·τ·S²** (S = singles rate), **front-loaded**
(∝ activity²), a few % of trues for τ ~ few ns. Method: time-tag each single from its parent
annihilation's time, sort by time, pair **cross-annihilation** singles within τ → `truth=2`.

## The blocker + the path

Randoms need a **real time** per single, from its **parent isotope's activity curve**. Today
times are dummy and the source is a geometric phantom with **no isotope** — real per-isotope
timing depends on Step 1 (the scenario source). **Path B (locked):** build the whole randoms
machinery now on a **toy activity model** (one isotope, single-exponential over an acquisition
window), wire real per-isotope activity later. Validation is the `2τS²` scaling, which the toy
model exercises fully.

## Locked decisions

- **Times are a downstream, deterministic, per-event quantity — NOT in the transport-only singles
  stack.** A pure function `event_time(model, ev, time_seed) -> Float64` (own seed, cheap
  per-event RNG e.g. `Xoshiro(hash(time_seed, ev))`), keyed by `event_number`. So: independent of
  the MT transport's chunking (alignment is by event id), reproducible over `(model, time_seed)`,
  O(1) memory (no 800 MB times-vector at 10⁸), recomputable identically in any pass.
- **A new single-purpose script** `scripts/build_randoms_from_singles.jl` (cross-annihilation
  pairing) — separate from `build_coincidences_from_singles.jl` (same-annihilation). Different
  algorithm (global time-sort vs event-streaming).
- **Build single-threaded first**, benchmark on the 10⁷ file, MT only what the numbers justify
  (see Benchmark below). Both consumers are currently single-threaded.
- **The complete measurement = trues+scatters ∪ randoms, time-ordered**, in one `lors.h5` with
  `has_randoms=true`.

## Pieces (ordered; ST first)

1. **`src/activity.jl`** — the toy activity model + `event_time`.
   - `event_time(model, ev, seed)`: sample from a truncated exponential `A₀e^{-λt}` over
     `[0, t_acq]` (λ = ln2/half_life), or flat. The singles rate `S = N_singles/t_acq` sets the
     randoms fraction → `t_acq` is the knob.
   - Real (later): per-isotope curves weighted by the scenario budget Nⱼ; `event_time` then reads
     the single's isotope tag (added by the Step-1 scenario source).

2. **`src/detector.jl`** — `smear_time(t, σ_t, rng)` (sibling of `smear_energy`);
   `[detector].sigma_t_ns`. The window-resolution blur on each single's time.

3. **`scripts/build_randoms_from_singles.jl`** (single-threaded first) — the randoms pass:
   - read singles, keep **good** ones (contained-one + pass the energy selection — reuse the
     shared `Response` / `pass_energy` / `fill_singles!` from `src/coincidences.jl`);
   - assign each its `event_time(ev)` (+ `σ_t` smear);
   - **sort by time**; scan the τ-window, pairing **cross-event** singles (same-event pairs are the
     trues — skip them: both photons of one annihilation share `t_ev`);
   - emit random LORs (`truth=2`) via `CoincidenceWriter` → `prod/<tag>/randoms_{mode}.h5`.
   - **Multiples policy** (design decision at impl): standard sorter rejects a window with >2
     singles; first cut may pair cross-event pairwise — note whichever is chosen.

4. **Real times into the trues/scatters LORs.** Extend `build_coincidences_from_singles.jl` to
   stamp each LOR's `t1,t2` with `event_time(ev)` (+ σ_t) instead of `0.0`. Cheap — it has `ev`.

5. **Merge → the complete `lors.h5`.** Combine trues+scatters + randoms, **sort by coincidence
   time**, write one file with `has_randoms=true`. A `merge_lors.jl`, or `build_randoms` appends.
   (At 10⁷ an in-memory concat+sort is fine; at 10⁸ this is another big sort — see MT below.)

## The pairing algorithm

- **Scale:** 10⁷ → ~13.7 M singles. First cut: **in-memory global sort** of the good singles'
  times (hundreds of MB — fine). For 10⁸: **τ-bin streaming** (bin into τ-width time bins, pair
  within/adjacent bins) — streaming-friendly and the natural MT partition.
- **Cross-event only** (same-event = trues). **No opposition filter** (randoms are accidental,
  any two blocks).

## Benchmark (on the current 10⁷ file) → the MT decision

The two consumers parallelize differently, so benchmark each:
1. **`build_coincidences_from_singles` (~5 s ST at 10⁷):** split into HDF5 read+decompress vs
   select+smear vs write (probe: time a read-only column load vs the full run). If
   **decompress/I/O-bound** → MT buys little (attack via lighter compression / chunked-parallel
   HDF5, not `@threads`); if **CPU-bound** → MT scales by **event-range chunks** (events
   independent — reuse the `simulate_source_mt` per-chunk-RNG + part-file ordered-merge pattern;
   note the smear RNG goes per-chunk → reproducible over `(seed, nchunks)`, not bit-identical to ST).
2. **`build_randoms` sort+pair:** prototype ST, measure sort + window scan, extrapolate to 10⁸.
   Sorting is CPU-bound → the likeliest MT win; parallelize by **τ-bins** (chunk by time) or a
   parallel sort.
3. **Deciding ratio:** each step's *ST 10⁸* time vs the transport's 10⁸ time (~30–60 s). Under
   ~10 s → not worth MT; sort-heavy randoms probably is.

## Validation

- **Rate** ≈ analytic `2τS²` (vary τ → linear; vary S via `t_acq` → quadratic).
- **Front-loaded** — randoms-vs-time ∝ activity² (falls faster than trues ∝ activity).
- **Magnitude** — a few % of trues for τ ~ few ns.
- Extend `scripts/tests/check_lors.jl` to report the three-way split (true/scatter/**random**) and
  expect `has_randoms=true` on the merged file.

## Config additions

```toml
[timing]   t_acq_s, half_life_s, shape (exp|flat), time_seed, tau_ns
[detector] sigma_t_ns
```

## Deferred

- **Real per-isotope activity** — needs Step 1 (scenario source: `emitters.csv` isotope tags +
  budget). The toy model unblocks everything else now; swapping it in later only changes
  `event_time` + adds an isotope column to the singles.
- **MT** for whichever consumer the benchmark justifies.
- A Julia reader/reducer for the production LOR HDF5 (analysis side).
