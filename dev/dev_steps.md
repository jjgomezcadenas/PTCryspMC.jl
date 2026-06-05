# Development steps

A running log of the simulation build. One entry per step: goal, what was done,
how to run it, status.

---

## Step 1 — phantom + 511 keV photon transport → photon stack

**Goal.** Define the phantom (a cylinder) from JSON, propagate 511 keV photons
through it, and write a CSV with the photon stack (interactions) per event.

**What was done.**
- `geometry/phantom.json` — the phantom: a water cylinder, radius 8 cm, half-length
  8 cm (the standard brain Ø16 × 16 cm), axis along z, centred at the origin.
- `data/materials.json` + `data/xcom_water.csv` — Water (density 1.0 g/cm³) with the
  NIST XCOM table (mass attenuation, cm²/g). `xcom_water.csv` is a NIST XCOM dump for
  H₂O including an exact 511 keV row.
- `src/` modules, adapted from LXeMC, trimmed to photon-only:
  - `nist_data.jl` — XCOM loader + log-log interpolation.
  - `materials.jl` — `Material`; `sigma_macro(mat, E)` returns the macroscopic cross
    sections Σ = (μ/ρ)·density [cm⁻¹] for incoherent/pair/photoelectric.
  - `geometry.jl` — `Cylinder`, `is_inside`, `distance_to_exit`, `load_phantom`.
  - `sampling.jl` — `sample_distance`, `sample_process`, `sample_compton`
    (Klein–Nishina), `rotate_to_global`.
  - `transport.jl` — `propagate_photon`: step through the cylinder; Compton (deposit
    recoil, continue) or photoelectric (deposit all, stop); record an `:escape` row at
    the exit point.
- `scripts/propagate_gammas_in_phantom.jl` — pencil source: 511 keV photons enter at
  the centre of the −z face along +z; writes `output/phantom_stack.csv`.

**Physics notes.**
- Only photons are followed; the recoil electron deposits locally (sub-mm range).
- Pencil along the axis: unscattered photons exit straight; scattered ones fan out,
  so the exit opening angle from Compton is recoverable from the stack.
- `mfp` in water at 511 keV ≈ 10.4 cm (μ ≈ 0.096 cm⁻¹), so an 8 cm half-length
  phantom is a fraction of a mean free path — most photons cross with 0–1 scatters.

**Output schema** (`output/phantom_stack.csv`, one row per interaction):
`event_number, step, x_mm, y_mm, z_mm, e_in_keV, e_dep_keV, process`
(`event_number` repeats across an event's rows; `process` ∈ compton / photoelectric /
below_cut / escape). The `escape` row carries the exit point (e_dep = 0).

**Run.**
```
julia --project=. scripts/propagate_gammas_in_phantom.jl --nevents 10000
```

**Status: done & validated.** 10⁴ pencil photons, water Ø16×16 cm.
- `mfp@511 keV = 10.44 cm` (μ ≈ 0.096 cm⁻¹) — correct for water.
- Unscattered fraction **0.215** vs analytic exp(−16/10.44) = **0.216** — the
  transport reproduces Beer–Lambert attenuation through the geometry.
- Each event terminates once (escape / photoelectric); scattered photons exit
  through the far cap, the side, and the entrance cap (backscatter), so the Compton
  opening angle is recoverable from the stack. Mean energy deposited ≈ 210 keV/event.
- `Pkg.test` passes; the module precompiles.

Next (step 2): the detector ring (block/wheel `CylShell`) and propagation of the
photon pair from an annihilation point to the ring.
