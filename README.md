# PTCryspMC.jl

Proton-therapy CRYSP simulations: a fast, photon-only Monte Carlo of a PET scanner's response
(in Julia).

The code produces, for a given scanner, the list of coincidences (lines of response) it would
record from positron annihilations, with each line tagged **true**, **scatter**, or **random**.
That list is the input to a separate range-precision analysis.

## Two modes

A shared physics engine (geometry, photon transport, detector response) is driven in two modes:

- **Phantom-based** — a geometric phantom (a sphere or cylinder of water) filled with a uniform
  source. The validated, currently-running mode; used to develop and check the engine, and a
  detector-study tool in its own right.
- **Proton-beam-based** — the eventual target: the positron activity a proton dose leaves
  (β⁺ emitters from nuclear fragmentation), drawn from a scenario produced upstream by `ptcryspg4`.
  Outlined; not yet built.

## Documentation

- `docs/PTCryspMC_phys.tex` — the engine: photon physics, geometry, transport, detector response.
- `docs/PTCryspMC_app.tex` — the application: the phantom-based and proton-beam-based modes.
- `CLAUDE.md` — orientation, decisions, and code layout. `dev/status.md` — current status and the
  deferred-work register.

The production chain (phantom mode): `simulate_source_mt.jl` → `singles.h5` →
`build_true_coincidences_from_singles.jl` + `build_randoms_from_singles.jl` → `reco_lors.jl` →
`lors_det.h5` (the list-mode deliverable).
