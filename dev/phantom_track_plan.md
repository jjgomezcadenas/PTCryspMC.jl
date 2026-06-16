# Plan — analytic-phantom → LOR-list track (Source, acollinearity, smearing, plotter)

Self-contained plan for the controlled-phantom validation track. Status at time of
writing: Step 3a (the multi-volume navigator) and a truth-level coincidence builder
(`scripts/build_coincidences.jl`) are done and committed. This track adds a proper
emission `Source`, a `Sphere` solid, detector smearing + the energy window, and a
coincidence plotter, so we can generate **list-mode coincidence files** from known
phantoms. No reconstruction here — the deliverable is the LM coincidence file.

## Goal

A controlled **analytic-phantom → list-mode coincidence file** chain: a phantom volume
(`solid + material`) that is *both* the source (uniform activity) *and* the attenuator,
emitted as back-to-back pairs **with acollinearity**, navigated to the ring, run through
**detector smearing + an energy window**, and written as a list-mode coincidence list (one
record per coincidence = one LOR), with a **diagnostics plotter**. Phantoms now:
**uniform cylinder + sphere** (Jaszczak/Derenzo later, on request).

## Locked decisions (from review)

- **Source = attenuator** (no decoupling): the phantom is a `solid + material` filled with
  uniform activity. The same volume both emits and attenuates; the `Source` samples
  uniformly inside `geom.phantom`'s **solid** (shape), independent of material.
- **Material is a free parameter** of the phantom. Filling the solid with **Vacuum/Air
  (Σ=0)** gives the **non-attenuated reference** — pairs emitted, nothing attenuated or
  scattered (LORs straight, every coincidence tagged `true`), recovering the pure source
  geometry modulo acollinearity + detector resolution. **Water** gives realistic
  attenuation + scatter. Comparing the two isolates the phantom-scatter contribution.
  (The navigator already handles a Σ=0 leaf: `propagate_photon` returns a single straight
  `:escape` segment, no interaction — no code change.)
- **Acollinearity** is on (~0.5° FWHM): the two photons are not exactly 180° apart.
- **Plotter** and **smearing** (Step 4) are both in scope.
- **True emission point** `(x0, y0, z0)` is carried in the stack and coincidence CSV — used
  only to **validate the source** (histogram to confirm the phantom is uniformly filled),
  not to feed any reconstruction.
- **No reconstruction at all.** The track's deliverable is the **list-mode coincidence
  file** (one record per coincidence = one LOR), which `build_coincidences.jl` already
  produces. Inverting LORs into an image (range precision, detector comparison) is the
  separate, deferred downstream analysis (CLAUDE.md: possibly its own repo).
- **Phantoms now:** uniform **cylinder** and **sphere** only. Jaszczak *and* Derenzo come
  later, on request (a composite `Source` of weighted sub-regions).

## Phase A — `Source` + `Sphere` solid (emission)

- `src/geometry.jl`: add **`Sphere(radius_cm) <: Solid`** with `is_inside`,
  `distance_to_entry`, `distance_to_exit`, `volume` (ray–sphere is a simple quadratic,
  `‖p + t d‖ = R`), and a `"sphere"` branch in `load_solid`; export it. The navigator
  needs **no change** — a sphere phantom is just another leaf absorber, reached through
  the existing `PhysicalVolume` delegation.
- `src/source.jl` (new):
  - **`UniformVolumeSource(pv::PhysicalVolume)`** — samples annihilation points uniformly
    inside the phantom's solid, via `sample_point_in(solid, rng)` for `Cylinder` (have
    the math), `Sphere` (`r = R·cbrt(u)`, isotropic), `Box` (uniform per axis), shifted
    by the volume's placement.
  - **`emit_pair(source, rng; acol_fwhm_deg=0.5)`** → `(pos, dir1, dir2)`: draw `pos` from
    the source, an isotropic `dir1`, and `dir2 = −dir1` tilted by a small Gaussian polar
    deviation (σ = FWHM/2.355, random azimuth) about `−dir1` — the acollinearity. Moves
    emission out of the driver script (the `Source` CLAUDE.md flagged).
- Exports: `Sphere`, `Source`/`UniformVolumeSource`, `emit_pair`, `sample_point_in`.
- Tests: sphere ray distances (entry/exit/inside, axial & off-axis, miss); uniform
  sampling stays inside the solid and fills it (rough density check); acollinearity angle
  distribution ≈ 0.5° FWHM (mean/σ within tolerance over many draws).

## Phase B — phantom geometry configs

- Add `geometry/geometry_sphere.json` (phantom = **sphere**); keep the existing
  `geometry.json` (**cylinder**). Choosing a phantom shape = choosing the geometry file
  (`--geometry`). Sphere radius defaults to **R = 8 cm** (matches the cylinder radius, for
  a direct comparison) — just a JSON value, editable anytime.
- The phantom **material** is a separate axis: a **`--phantom-material`** override on the
  driver (default = the JSON value, like the existing scanner `--material`) selects
  Vacuum/Air (non-attenuated reference) vs Water (realistic) without new geometry files.

## Phase C — driver

- Generalize `scripts/simulate_phantom.jl` to emit via
  `UniformVolumeSource(geom.phantom)` + acollinearity, so `--source phantom` works for a
  cylinder *or* a sphere phantom, writing the same tagged stack.
- Add the **`--phantom-material`** override (default = JSON) so the same run can be done
  with a Vacuum/Air phantom (non-attenuated reference) or a Water phantom (realistic).
- Optionally add the **true emission point** `(x0, y0, z0)` to the stack / coincidence CSV
  — carried so the plotter can confirm the source is uniformly filled (source
  validation), not for reconstruction.

## Phase D — smearing + energy window in `build_coincidences.jl` (Step 4)

Extend the streaming builder with optional detector response. **Truth mode (no flags) is
unchanged** (the current behaviour). Detector mode adds:
- `--sigma-xyz` (mm, default **1.7**, CRYSP) → 3-D Gaussian smear of each formed hit
  position (includes DOI, since the hit is the first-interaction 3-D point);
- `--eres a` (FWHM at 511 keV, **0.05** CsI / 0.07 CsI(Tl) / 0.10 BGO) → smear each hit
  energy by FWHM(E) = a·√(511 keV / E);
- `--window` (half-width in FWHM units, default **2**) → keep the pair iff both smeared
  energies fall in 511 ± window·FWHM(511 keV);
- `--seed` for the smear RNG.

This converts the loose topological list into the realistic one: low-energy scatters drop
out at the window, near-full scatters survive tagged `scatter`. **Time smear σ_t is
deferred** — there is no real `t` until the randoms pass, so `t` stays the dummy 0.0.

## Phase E — coincidence plotter (Python, `py/plot_coincidences.py`)

Reads a coincidence CSV — **diagnostics only, no reconstruction**. Panels: energy spectra
`e1`/`e2` with the window overlaid; the true/scatter split; the φ–z endpoint map; a radial
profile of the hits; and the **true source distribution** from `(x0,y0,z0)` (x–y and x–z
density) to confirm the phantom is uniformly filled. (Plotting stays Python per the updated
CLAUDE.md tech-stack.)

## Phase F — run + check

Run the full chain for the cylinder and sphere phantoms (CsI and BGO), both Vacuum
(non-attenuated reference) and Water: confirm the source `(x0,y0,z0)` fills the phantom
uniformly, the Vacuum case is all-`true`, and the true/scatter split behaves sensibly as
the energy window tightens. The deliverable is the list-mode coincidence file.

## Files to touch

- `src/geometry.jl` — `Sphere` solid + `load_solid` branch + export.
- `src/source.jl` (new) — `Source` types, `sample_point_in`, `emit_pair` (acollinearity).
- `src/PTCryspMC.jl` — include + exports.
- `test/runtests.jl` — sphere, uniform sampling, acollinearity testsets.
- `geometry/geometry_sphere.json` (new).
- `scripts/simulate_phantom.jl` — emit via the `Source`; optional `(x0,y0,z0)`.
- `scripts/build_coincidences.jl` — smearing + energy window (Step 4 flags).
- `py/plot_coincidences.py` (new) — the coincidence plotter.

## Choices — resolved

1. **Sphere radius** — default **R = 8 cm** in `geometry_sphere.json` (just a JSON value).
2. **Phantoms** — uniform **cylinder + sphere** now; **Jaszczak and Derenzo both later**,
   on request.
3. **True emission point `(x0,y0,z0)`** — **yes**, carried for source validation only.
4. **Reconstruction** — **none**. The deliverable is the list-mode coincidence file;
   inverting LORs to an image is the separate, deferred downstream.

## Suggested order

Phase A (Source + Sphere) first — it is foundational and testable in isolation. Then B/C
(geometry + driver), D (smearing + window), E (plotter), F (validate). Jaszczak/Derenzo is
a later addition to the `Source` (a composite of weighted sub-regions).
