# Plan — Step 5: randoms (and the timing it needs)

Self-contained plan so this can be picked up after a context clear. Status: Steps 1–4 of
the pipeline are effectively done via the analytic-phantom → LOR track (A–F, see
`dev/phantom_track_plan.md`) — transport, navigation, hit formation, smearing, the energy
selection, the list-mode coincidence file (true/scatter), config-driven and parallel.
**Times are still a dummy `0.0`** everywhere. Randoms is the next major piece.

## What a random is

Two photons from **different** annihilations that happen to arrive within the coincidence
window τ (a few ns) → a false coincidence with a wrong LOR. Rate ≈ `2·τ·S²` (S = singles
rate), **front-loaded** (largest at scan start, falls as the activity²), at most a few % of
the trues for τ ~ a few ns. Method: a **separate pass over the singles, no re-transport**
(`docs/pet_simulation.tex`): time-tag each single from its isotope's activity, sort in
time, pair cross-annihilation singles within τ, tag them `random`.

## The blocker: randoms needs real TIMES, which need ISOTOPES

Each single needs an absolute time drawn from its **parent isotope's activity curve**
(e.g. ¹⁵O, 122 s half-life; front-loaded). Today:
- times are dummy `0.0` (no timing in the sim);
- the source is a **geometric phantom** (`UniformVolumeSource`) with **no isotope** — so
  there is no activity curve to sample from.

So real randoms depends on **Step 1 — the scenario source** (read `ptcrysp-scenarios`:
`emitters.csv` = annihilation points + isotope per emitter, `sampling_budget_*.csv` = N_j
per isotope, `run_meta.csv` = dose/timing, `isotopes.csv`). That tags each annihilation
with an isotope → its activity curve → its time.

**Two paths:**
- **Path A — scenario source first (Step 1), then randoms.** The "correct" order. Step 1 is
  its own chunk: read the scenario, draw N_j isotope-tagged annihilation points to the
  budget, emit. Then time-tag from per-isotope activity, then randoms.
- **Path B — simplified randoms first.** Build the randoms machinery (time-tag → sort →
  pair-within-τ) on the current phantom source with a **toy activity model** (one isotope,
  a single exponential `A(t)=A₀e^{-λt}` over the acquisition window, or a flat rate), to
  develop + test the pass; plug in the real scenario activity later. *(Recommended to start
  — decouples the randoms algorithm from the scenario I/O.)*

## What needs building

1. **A singles list.** Per *detected photon* (one row per single): position, energy, block,
   **time**, parent isotope, scatter flag. This is the Phase-G `--singles` reduced stack
   (`dev/phantom_track_plan.md` Phase G) **plus** an isotope tag and a sampled time. So
   Step 5 builds on the singles stack — do that part of Phase G alongside.
2. **Timing.**
   - **σ_t** (time resolution) → add `smear_time` to `src/detector.jl` (sibling of
     `smear_energy`/`smear_position`); `[detector].sigma_t_ns` in the config.
   - **Activity model** → sample each annihilation's time from its isotope's decay curve
     over the acquisition window (Path B: toy model; Path A: from the scenario).
3. **The randoms pass** (`scripts/build_randoms.jl`): read the singles, time-tag, and pair.
   - Pairing needs **time order** — *not* a single streaming pass like coincidences. At 10⁸
     singles a global sort is heavy; options: (a) in-memory sort if it fits, (b) **time-bin**
     into τ-width bins and pair within adjacent bins (streaming-friendly). Decide at impl.
   - For each cross-annihilation pair within |Δt| ≤ τ that passes the energy selection and
     forms a valid two-block LOR → emit a coincidence tagged `random`.
4. **Output.** Append randoms (`truth=random`) to the coincidence list, or a separate
   `randoms_*.csv`. The plotter's truth split then has three categories (true/scatter/random).

## Design questions to resolve when we start

1. **Path A or B** (scenario source first, or toy-timing randoms first)? — lean **B**.
2. **Singles list**: build it via the Phase-G `--singles` stack + isotope + time?
3. **Pairing at scale**: global sort vs τ-bin streaming?
4. **Output**: randoms appended to the coincidence file, or separate?
5. **Activity model** location: a small `src/activity.jl` (per-isotope curves, half-lives)?

## Adjacent deferred items (related)

- **σ_t time smearing** — only meaningful once times exist (this step).
- **Phase G remainder** — `--singles` reduced stack + HDF5; the singles list is the natural
  place these meet.
- **Step 1 — scenario source** — the prerequisite for real (isotope-based) timing; itself a
  piece of work (ptcrysp-scenarios I/O), independent of the randoms *algorithm*.

## Quick state reference (for resume)

Pipeline today (all committed, `main`): `runs/<tag>.toml` → `simulate_phantom.jl`
(transport → `output/<tag>/stack.csv`, with dummy times) → `build_coincidences.jl`
(streaming select + smear + energy cut → `coincidences_{truth,det}.csv`, `truth ∈
{true,scatter}`) → `plot_coincidences.py` / `plot_matrix.py`. Launchers:
`scripts/run_matrix.sh` (data, parallel), `scripts/plot_all.sh` (plots, parallel).
`Pkg.test`: 202. Detector physics in `src/detector.jl`; emission in `src/source.jl`;
config in `src/config.jl`.
