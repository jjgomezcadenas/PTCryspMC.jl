# Project status — PTCryspMC.jl

A one-page snapshot of where the simulation stands, plus the **deferred-work register**.
Companion docs: the *method* is in `docs/pet_simulation.tex` and `docs/overview.tex`; the
decisions + code layout in `CLAUDE.md`; the full chronological build log in `dev/dev_steps.md`.

_Last updated: 2026-06-24._

## What it does

A fast, photon-only Monte Carlo: read a frozen scenario (positron-annihilation points left by a
proton field) + a detector description → write the list-mode coincidence/LOR file a PET scanner
would record. It never runs proton transport.

## The pipeline (built & validated)

```
runs/<tag>.toml
   │
simulate_source_mt.jl  (multi-threaded, alloc-free) ──► prod/<tag>/singles.h5   (+ t_rel)
   ├─ build_true_coincidences_from_singles.jl ──► lors_truth.h5   (true + scatter; t1,t2,dt)
   ├─ build_randoms_from_singles.jl            ──► randoms.h5      (truth = random)
   └─ reco_lors.jl  (smear + energy + DT cut + flag) ──► lors_det.h5   (the list-mode deliverable)
```

- **Timing model.** Each single carries `t_rel` = TOF + scintillation first-photon jitter
  (`jitter = −ln u / (N_det·r0)`, `N_det = yield·E·pde`, `r0 = Σ wₖ/τₖ`), stamped **once** at
  generation, **relative to the decay** so it stays small (Float32). Absolute time =
  `event_time(ev)·1e9 + t_rel`, reconstructed only where randoms need a common clock.
- **Coincidence window τ = 3 ns** (CsI and BGO), read off the DT study (`examine_dt.jl` +
  `py/plot_dt.py`); `compare_crystal_timing.jl` explains why τ is crystal-independent.
- **Reco lower energy cut** `reco_emin_keV = 450` (photopeak region; the spectrum studies keep
  `emin_keV = 300` to see the Compton shoulder). No upper cut (it's an analysis-time choice).
- **Validated (10⁷):** `Pkg.test` **812**; randoms match the analytic `2τS²` (CsI 1248 vs 1291,
  ratio 0.97); reco acceptance CsI 8.98% / BGO 23.75%; corrected residual median|dt| (CsI 0.059 /
  BGO 0.200 ns) matches the analytic single-photon jitter.
- **Scale.** `scripts/tests/bench_chain.jl`: the full chain at 10⁸ runs **~3 min serial, ~14 GB
  peak** (`build_randoms`, the only N-scaling stage) — fits 48 GB, no rework needed.

## Detector configs

CsI (5% FWHM @ 511 keV) and cryogenic BGO (10%) — sphere-water scenarios, validated. CsI(Tl) (7%)
and the pixelated detectors (LYSO, standard BGO) are still to add.

## Documentation

- `docs/pet_simulation.tex` — the method (detailed).
- `docs/overview.tex` — the project overview (physics + the geometry/source/detector choices).
- Docstrings + a Literate/Documenter web doc-site — **planned (the next doc task), not built.**

## Deferred work & known nits

From the multi-agent code review. **None are correctness bugs** — all are robustness/clarity.
Each numbered item has a natural trigger that makes it cheap to fold in *then* rather than now;
the corresponding code site carries a short `Deferred:` marker pointing here.

### Deferred fixes (specific trigger)

| # | Where | What | Trigger |
|---|-------|------|---------|
| 1 | `src/timing.jl` `first_photon_jitter` | `-log(rand)` → `+Inf` if `rand()==0` (p ≈ 2⁻⁵³); fix `-log(1-rand)`. | Changes every jitter value → **forces a singles regeneration**; fold in next time the singles are regenerated for another reason. |
| 2 | `src/geometry.jl` `Scanner.volume` | the field named `volume` (a `PhysicalVolume`) collides with the generic `volume()` (cm³); rename → `pvol`/`shell`. | A small rename sweep (`geometry.jl` + `simulate_source_mt.jl`); do it when next touching the Scanner API. |
| 3 | `src/nist_data.jl` `load_xcom` | the parser assumes 8 clean numeric columns + digit-first data rows; XCOM **K-edge label rows** can shift columns, and empty input throws obscurely. | Adding the **CsI/BGO/LYSO XCOM tables** (which have in-range K-edges) — harden with the real tables to test against. |

### Minor nits (cosmetic / robustness, no trigger)

- `block_index` docstring (geometry): φ ∈ [0,2π) never equals 2π (the clamp is FP-defensive); `z < −H` silently clamps to `iz=0`.
- `_prepare_xcom_energy` (nist_data): the K-edge nudge assumes the below-edge row precedes the above-edge row (true for standard XCOM) — note it when crystals arrive.
- `_leaf_reduce` (navigator): the inner Compton loop has no explicit step backstop; it terminates because energy strictly decreases each Compton to the cut — worth a one-line comment.
- `position_cm` (geometry loader): a JSON `position_cm` with ≠3 entries throws a confusing `TypeAssertionError`; add a `length==3` check.
- HDF5 chunk size: a fixed 2²⁰-row chunk regardless of run size (oversized/zero-padded for tiny runs); the `total==0` empty-run path is fine but undocumented.
- `build_randoms` Float32 edge: `t2 = Float32(t_rel+Δ)`, then reco recomputes `|t1−t2|`; Float32 rounding can nudge a `Δ ≈ τ` random just across the cut. Sub-permille.

### Deferred by scope (not from the review)

- Real per-isotope activity (Step 1 scenario source: isotope tags + the sampling budget). Swapping it in changes only `event_time` + adds an isotope column to the singles.
- Refine the per-crystal **PDE** (0.45 placeholder for CsI and BGO; should differ by emission colour).
- Cylinder/vacuum configs need their own `examine_dt` (different scanner length → different TOF tail) before fixing their `[timing].tau_ns`.
- Threshold / CFD timing (the "first photon" model is the leading-edge idealization).
- Pixelated detectors report a fixed crystal, not a continuous position. Placement **rotation** transform (only when a volume needs it).

## Next

1. **Docstrings + Literate/Documenter doc-site** (the remaining documentation task).
2. Crystal XCOM tables + the R5 `load_xcom` hardening; the CsI(Tl) config.
3. Downstream, separate & deferred: range precision, detector comparison, reconstruction.
