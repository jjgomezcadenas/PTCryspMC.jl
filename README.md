# PTCryspMC.jl

Proton-therapy CRYSP simulations: a fast, photon-only Monte Carlo of a PET scanner's response
(in Julia).

The code produces, for a given scanner, the list of coincidences (lines of response) it would
record from positron annihilations, with each line tagged **true**, **scatter**, or **random**.
That list is the input to a separate range-precision analysis.

## Two modes

A shared physics engine (geometry, photon transport, detector response) is driven in two modes,
selected by `[source].mode` in the run config:

- **Clinical** (`mode = "clinic"`) — a radiotracer distribution at a known activity: one or more
  regions (a named geometry volume, or an inline shape), each with a concentration or a total
  activity; an isotope (F-18 default) and an acquisition window derive the number of annihilations
  from the decay curve. Structured phantoms (e.g. NEMA-IQ) need no geometry change. **Built and
  validated** — the currently-running mode; the count-driven uniform phantom survives as its
  special case. Detectors so far: pure CsI and cryogenic BGO.
- **Proton Activity** (`mode = "api"`) — the positron activity a proton dose leaves behind (β⁺
  emitters from nuclear fragmentation), read from a frozen `ptcryspg4` scenario with its per-isotope
  decay budgets. **Built and validated.** Runs export to the `PtCryspProds/` products tree
  (`scripts/run/publish_prod.jl`) for the downstream reconstruction / range-precision repo.

## Running

- **Production chain** — per *named* config (it never runs unnamed):
  `scripts/run/run_prod.sh <config> [config ...]` → `prod/<tag>/`. The chain is
  `simulate_source_mt.jl` → `singles.h5` → `build_true_coincidences_from_singles.jl` +
  `build_randoms_from_singles.jl` → `reco_lors.jl` → `lors_det.h5` (the list-mode deliverable),
  plus a control-plot PNG.
- **Dev chain** (full interaction stacks, CSV): `scripts/run/run_matrix.sh` → `output/<tag>/`.
- **Tests**: `julia --project=. -e 'using Pkg; Pkg.test()'`.

## Documentation

- `docs/PTCryspMC_phys.tex` — the engine: photon physics, geometry, transport, detector response.
- `docs/PTCryspMC_app.tex` — the application: the Clinical and Proton Activity (API) modes.
- `docs/SCHEMA.md` — the output-file schema (generated from the code; a test keeps it in sync).
- `CLAUDE.md` — orientation, decisions, and code layout. `dev/status.md` — current status and the
  deferred-work register.

Build the manual PDFs with `python3 py/build_latex.py` (compiles `docs/*.tex` and clears the LaTeX
aux clutter; figures are regenerated from `py/fig_*.py` with `--figures`).
