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

**Done.** Six `runs/*.toml` (Vacuum for CsI only), 100k events each, run in parallel via
`scripts/run_matrix.sh` (data) then `scripts/plot_all.sh` (plots); ~7 s wall for the data,
~1 s for the plots. The energy cut is a 300 keV **minimum** (`[detector].emin_keV`, window
off) so the Compton shoulder/edge shows. Validated: Vacuum = 100 % true; BGO ~2.6× the CsI
yield; LOR resolution ~1.7 mm (true) across all configs (acollinearity + σ_xyz, material-
independent). Per-config 9-panel summaries (`plot_coincidences.py`) + a cross-config
comparison (`plot_matrix.py` → `output/matrix_summary.png`). A near-pole `rotate_to_global`
bug (tuple vs `Vector{Float64}` for a photon emitted ≈ along z, ~2/100k events) was found
and fixed + regression-tested.

## Phase G — production I/O: singles stack, CSV/HDF5, TOML config

CSV + the full per-interaction stack is ideal for the testing phase and detailed plots
(good to ~10⁵ events). Production (10⁸ decays) needs a compact stack and a binary format.

### The `--singles` (reduced) stack

The full stack is one row per *interaction* (~3–6 per photon) and is the size bottleneck
(order 10⁹ rows at 10⁸ decays). The **singles** stack is **one row per *detected photon***
— the formed hit — which is exactly the reduction `build_coincidences` already computes
internally (`GammaAcc`). It is the *singles list* the method calls for (randoms + the hit
truth), minus the per-interaction detail (which only the CNN R&D needs, run at small scale
on the full stack).

Singles schema (one row per photon that deposited in the ring; a miss writes nothing):
```
event_number, gamma, x_mm, y_mm, z_mm, e_keV, iz, iphi, nblocks, phantom_scatter, x0_mm, y0_mm, z0_mm
```
`x,y,z` = first scanner interaction (the LOR point); `e_keV` = summed crystal energy
(truth, unsmeared — smearing stays in `build_coincidences`); `nblocks` = distinct blocks
touched (1 = contained, >1 = overspill). ~5–10× fewer rows and narrower than the full stack.

`simulate_phantom.jl` writes the singles row directly when `--singles` is set (transport
is unchanged; only the I/O is reduced — the giant full stack is never materialised).
`photon_stats` is extended to also return the first scanner interaction `(x,y,z,iz,iphi)`.

### CSV / HDF5 (`--format csv|hdf5`)

Binary HDF5 for production: compact (compression), typed, partial/chunked reads (randoms
and the CNN pull just the events/columns they need), and self-describing — the run
parameters live in root **attributes** (the `_meta` companion, in the same file).
Streaming write = chunked extensible datasets (buffer N rows → append a chunk), so both
scripts keep O(1) memory. Adds `HDF5.jl` (Julia) + `h5py` (Python).

`build_coincidences` **auto-detects**: input format by extension (`.csv`/`.h5`), input
*kind* by header columns (full has `step`/`process`/`volume`; singles has `e_keV`/`nblocks`).
The selection + smearing (`finish_event!`) is **shared** — both paths just fill `GammaAcc`
(full: accumulate scanner deposits; singles: one row already is the hit). Its own LOR
output also takes `--format`.

### TOML config (the parameter source of truth) + run/output layout

A single pipeline config with sections, read by **both** scripts (each uses its relevant
sections), replaces the long CLI and **doubles as the run metadata**:
```toml
[geometry]  file, phantom_material
[source]    kind (point|phantom), acol_fwhm_deg, energy_keV, nevents
[transport] crystal_material, cutoff_keV, seed
[detector]  sigma_xyz_mm, eres, window_fwhm, seed
[output]    dir, tag, format (csv|hdf5), singles (bool)
```

**Directory layout.**
- `runs/` (**tracked** in git) — one config per condition, named `<shape>_<phantom_mat>_<crystal>.toml`,
  e.g. `runs/sphere_water_csi.toml`, `runs/sphere_vacuum_csi.toml`. These are the
  reproducible recipes (the provenance of every result), so they are versioned.
- `output/<tag>/` (**gitignored**) — per-run subdirectory, self-describing:
  `stack.csv`, `coincidences_{truth|det}.csv`, the plot `.png`, and a **copy of the
  config** (`config.toml`). With many conditions this keeps results from colliding.

**The `tag`** = `[output].tag` if set, **else the config filename** (without `.toml`). So
`runs/sphere_water_csi.toml` → tag `sphere_water_csi` → everything lands in
`output/sphere_water_csi/`. One token names the config, the output dir, and the metadata.

**Invocation.** `julia … simulate_phantom.jl --config runs/sphere_water_csi.toml` and
`build_coincidences.jl --config runs/sphere_water_csi.toml` (derives the stack from the
tag). TOML is the source of truth; a thin `--key value` override (e.g. `--nevents`,
`--seed`) is allowed for sweeps. Config reading uses the Julia `TOML` stdlib
(`src/config.jl`: `read_config`, `run_tag`, `cfg_get`); paths resolve against the repo
root so runs work from anywhere.

### Matrix
| scale | stack | format |
|---|---|---|
| dev / plots (≤10⁵) | full | csv |
| medium / CNN sample | full | hdf5 |
| production (10⁸) | singles | hdf5 |

Largely independent of E/F: the plotter reads whichever format/kind it is given, so the
dev chain can stay CSV+full while G is built. Do G when the scale demands it.

## Files to touch

- `src/geometry.jl` — `Sphere` solid + `load_solid` branch + export.
- `src/source.jl` (new) — `Source` types, `sample_point_in`, `emit_pair` (acollinearity).
- `src/PTCryspMC.jl` — include + exports.
- `test/runtests.jl` — sphere, uniform sampling, acollinearity testsets.
- `geometry/geometry_sphere.json` (new).
- `scripts/simulate_phantom.jl` — emit via the `Source`; `(x0,y0,z0)`; `--singles`, format, TOML.
- `scripts/build_coincidences.jl` — smearing + window; auto-detect format/kind; TOML.
- `py/plot_coincidences.py` (new) — the coincidence plotter (reads CSV or HDF5).
- **Phase G:** `Project.toml` (+ `HDF5.jl`); a shared CSV/HDF5 row-writer + TOML config
  reader (e.g. `src/io.jl`); a `run.toml` example.

## Choices — resolved

1. **Sphere radius** — default **R = 8 cm** in `geometry_sphere.json` (just a JSON value).
2. **Phantoms** — uniform **cylinder + sphere** now; **Jaszczak and Derenzo both later**,
   on request.
3. **True emission point `(x0,y0,z0)`** — **yes**, carried for source validation only.
4. **Reconstruction** — **none**. The deliverable is the list-mode coincidence file;
   inverting LORs to an image is the separate, deferred downstream.

## Suggested order

Phase A (Source + Sphere) first — foundational and testable in isolation. Then B/C
(geometry + driver), D (smearing + window), E (plotter), F (validate). **A–D are done and
committed.** Phase G (production I/O: `--singles`, CSV/HDF5, TOML) is largely independent
and can be built when the scale demands it; the dev chain stays CSV+full meanwhile.
Jaszczak/Derenzo are a later addition to the `Source` (a composite of weighted sub-regions).
