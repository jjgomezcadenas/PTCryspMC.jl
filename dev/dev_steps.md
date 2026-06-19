# Development log

The running record of the build: what the simulation is meant to do, what works
so far, and what comes next. The *method* lives in `docs/pet_simulation.tex`; this
file tracks the *implementation*.

---

## 1. The simulation

A fast, photon-only Monte Carlo that turns an upstream scenario (positron
annihilation points left by a proton field) into the coincidence list a given PET
scanner would record. It reads a scenario and a detector description and writes
`coincidences_<config>.csv`; it never runs proton transport.

The end-to-end pipeline:

1. **Source injection** — draw annihilation points per isotope to the measured
   budget; emit two back-to-back 511 keV photons each, with ~0.5° FWHM
   acollinearity.
2. **Transport** — follow only the photons through phantom (water) → air → crystal
   ring: free path from XCOM cross sections, Compton + photoelectric, electrons
   deposited locally. → a singles list + same-annihilation coincidences.
3. **Hit formation** — group deposits per detector block; position = first
   interaction point, energy = sum, time = earliest; smear by σ_xyz (incl. DOI),
   σ_t and FWHM(E) = a·√(511 keV / E); apply the energy window.
4. **Selection** — exactly two roughly-opposite blocks → a coincidence (true if
   neither photon scattered in the phantom, else scatter); anything else fails.
5. **Randoms** — a separate pass over the singles (no re-transport): time-tag each
   from its isotope's activity, pair cross-annihilation singles within the window τ.
6. **Detector configs** — the monolithic crystals: CsI (a = 5%), CsI(Tl) (7%),
   cryogenic BGO (10%), σ_xyz ~1.7 mm.

Resting under all of it: a **photon physics core** (XCOM cross sections; Compton
and photoelectric samplers; the step loop) and a **Geant4-style geometry**
(solids → logical volumes → physical volumes, placed in a world).

Downstream and deferred (separate, possibly its own repo): range precision,
detector comparison, reconstruction.

---

## 2. Achieved so far

Today the code loads a geometry and materials and **navigates** 511 keV photons
across the full geometry — phantom (water) → air → crystal ring, switching material
at each boundary — writing the per-interaction stack tagged with volume and detector
block. The transport (pipeline stage 2) is in place and validated against Beer–Lambert;
what remains is the scenario source (stage 1) and the coincidence/hit/randoms/detector
bookkeeping on top (stages 3–6).

### Foundations — geometry, materials & physics core, tooling

The reusable substrate the rest of the pipeline plugs into. (Built out while
reviewing the first transport pass; it is the base, not a detour.)

- **Geometry hierarchy** (`src/geometry.jl`), Geant4 semantics:
  - `Solid` (abstract) — pure shape in a local frame centred at the origin.
    `Cylinder <: Solid` implements `is_inside`, `distance_to_entry`,
    `distance_to_exit`, `volume`. The corner-case fix (inclusive `<=` boundary so a
    ray through the exact rim is never dropped), named tolerances (`SURFACE_EPS`,
    `PARALLEL_EPS`) and the shared allocation-free `_surface_crossings` scan live here.
  - `LogicalVolume{S<:Solid}` — `name` + `solid` + resolved `material`, with
    `volume` / `mass`. The reusable, placement-free level.
  - `PhysicalVolume{S<:Solid}` — a logical volume placed at a `position` in the
    world frame; the world→local transform (`_to_local`) is written once here and
    delegates to the solid (rotation slots into the same place when needed).
    Parametric `{S}` keeps every level concretely typed → the nesting is free at
    runtime (the `bench_geometry_levels.jl` probe: ~0.45 ns/call, 0 B, ratio ≈ 1.0).
  - **Full-geometry file.** `geometry/geometry.json` describes the world, one named
    section per component (`phantom` now, the detector ring later).
    `load_geometry(path, materials)` returns a `Geometry` container (so far the
    phantom physical volume) — the seed of the volume list a multi-volume navigator
    will walk. `load_solid` is the shape factory; unknown shapes and missing
    materials are rejected, not silently mis-loaded.
- **Materials & cross sections** (`src/nist_data.jl`, `src/materials.jl`):
  - `Material` + `sigma_macro(mat, E)` → the macroscopic cross sections
    Σ = (μ/ρ)·density [cm⁻¹], ordered `(Compton, photoelectric, pair)` (pair last,
    being zero below threshold).
  - `load_material(dir, name)` (single material — an unfinished entry can't break an
    unrelated run) and the batch `load_materials`, sharing one `_build_material` path.
  - `data/xcom_water.csv` — NIST XCOM for H₂O, 10 keV → 10 MeV, with an exact 511 keV
    row; the pair channel is real above ~1.5 MeV.
- **Physics core** (`src/sampling.jl`, `src/transport.jl`):
  - `sample_distance`, `sample_process` (photoelectric is the catch-all so a
    zero-width pair bucket can't absorb a rounding leftover), `sample_compton`
    (Klein–Nishina), `rotate_to_global`.
  - `propagate_photon(E0, pos, dir, pv, rng)` — steps through a physical volume's
    solid, reading its material; Compton (deposit recoil, continue), photoelectric
    (deposit all, stop), `:below_cut` and `:escape` records. Single-volume for now.
- **Tooling / control plots**:
  - `scripts/water_xsections.jl` → `output/water_xsections.csv` — the macroscopic
    cross sections at ~20 log-spaced energies (decoupled from the XCOM grid; refuses
    to extrapolate).
  - `py/plot_water_xsections.py` → `output/control_plots/water_xsections.png` — the
    log–log channel + total plot.
  - `scripts/bench_geometry_levels.jl` — the alloc-free / no-cost-per-level guard.

### First result — phantom + 511 keV transport → photon stack

The first running milestone: define the phantom from JSON, propagate 511 keV
photons through it, write the photon stack per event.

- `scripts/propagate_gammas_in_phantom.jl` — pencil source: 511 keV photons enter at
  the centre of the phantom's −z face along +z; writes `output/phantom_stack.csv`.
- **Output schema** (one row per interaction):
  `event_number, step, x_mm, y_mm, z_mm, e_in_keV, e_dep_keV, process`
  (`process` ∈ compton / photoelectric / below_cut / escape; the `escape` row carries
  the exit point with e_dep = 0).
- **Run:**
  ```
  julia --project=. scripts/propagate_gammas_in_phantom.jl --nevents 10000
  ```
- **Validated.** 10⁴ pencil photons, water Ø16×16 cm:
  - mfp@511 keV = 10.44 cm (μ ≈ 0.096 cm⁻¹) — correct for water.
  - Unscattered fraction **0.215** vs analytic exp(−16/10.44) = **0.216** — the
    transport reproduces Beer–Lambert through the geometry (and still holds after the
    geometry-hierarchy refactor).
  - Scattered photons exit through the far cap, the side, and the entrance cap
    (backscatter), so the Compton opening angle is recoverable from the stack; mean
    energy deposited ≈ 210 keV/event.

### Single crystal — CsI & BGO containment study

Characterising the crystal response (energy containment, interaction depth) on a
single crystal modelled as a box, before the full ring.

- **`Box` solid** (`src/geometry.jl`) — axis-aligned, analytic slab method for
  `is_inside` / `distance_to_entry` / `distance_to_exit` / `volume`; a `"box"`
  branch in `load_solid`.
- **Crystal materials** (`data/materials.json` + XCOM tables, 10 keV–10 MeV): CsI
  (4.51 g/cm³, attenuation length 2.44 cm @511 keV) and BGO (Bi₄Ge₃O₁₂,
  7.13 g/cm³, 1.10 cm), both with K-edges (I 33 keV, Bi 90.5 keV).
- **Cross-section tooling, now generic over material**: `material_xsections.jl`
  (was `water_xsections.jl`) + `py/plot_material_xsections.py` →
  `output/<material>_xsections.{csv,png}`.
- **`scripts/shoot_gammas_to_crystal.jl`** — 511 keV photons into a crystal box,
  fully parametric beam (`--beam-xy`, `--beam-opening` = total cone opening,
  `--material`, `--tag`) → `output/<material>_crystal_<tag>_stack.csv`.
  **`py/plot_crystal.py`** — a 3×3 panel (containment, photo/compton split, Etot,
  scatters, 1st/2nd Compton depth and separation, x–y interaction map).
- **Results** (48×48×37 mm, 511 keV; cone = 45° opening, pencil = on-axis):

  | beam | CsI contained | BGO contained |
  |---|---|---|
  | cone   | 48% | 84% |
  | pencil | 58% | 92% |

  BGO's short attenuation length and high-Z Bi give far higher containment and a
  larger photoelectric share (photo ≈ 48% of contained vs ≈ 33% for CsI).

### Detector geometry — air world + the scanner ring

The world model and the detector ring (Step 2 geometry).

- **Air world (mother volume)** — an `Air` material (no XCOM → zero cross section →
  straight-line propagation) and a `world` cylinder of Air (radius 60, half-length
  60 cm) in `geometry.json` that encloses every daughter. `Geometry` is now
  `world` + `phantom` + optional `scanner`; the radial layering is water [0,8] →
  air → crystal [38.7,42.4] cm.
- **`CylShell` solid + `Scanner`** (`src/geometry.jl`) — the hollow-cylinder ring,
  with the non-convex ray intersection by the interval method
  (S = (outer ∩ z-slab) ∖ bore, up to two pieces). The `(φ, z)` block/wheel grid:
  `block_index` / `block_id` (= iz·n_phi + iφ) / `nblocks`. The `scanner` daughter
  is `CRYSP_CSI_1M` (CsI: Ri=38.7, wall=3.7, H=51.2 cm, 48 φ × 20 z = 960 crystals).
- **`docs/navigation.tex`** documents the ray–cylinder / ray–shell distances and the
  (φ, z) partition (φ-pitch ≈ 50.6 mm, z-pitch 51.2 mm ≈ the 50 mm crystal).

### Unit test — back-to-back 511 keV pairs into the ring

A point source emitting back-to-back pairs isotropically into the ring — the first
source→ring path. The interior is air (no phantom scatter), so the navigation is
trivial: straight through air to the crystal (`distance_to_entry`), then transport.

- **`scripts/shoot_into_ring.jl`** — script-local `emit_pair`
  (no phantom/scenario), settable low-energy cutoff (default 10 keV = XCOM min).
  One CSV records both photon stacks tagged with the crystal `(iz, iphi)`:
  `event_number, gamma, step, x,y,z, e_in, e_dep, process, iz, iphi`.
- **`py/plot_back_to_back.py`** — a 3×3 panel (Edep, 3-D impacts, DOI of the 1st/2nd
  interaction, per-coincidence containment, contained-Edep, impact within a crystal
  face, overspill, φ–z hit map). **`scripts/bench_back_to_back.jl`** — ~0.5 µs/event.
- **Results** (point source, CRYSP1M ring):
  - ring acceptance **79.5%** of photons — geometric, material-independent
    (|cosθ| ≤ 38.7/√(38.7²+51.2²) = 0.798).
  - per-gamma containment / clean coincidence (both γ contained), CsI **0.54 / 0.29**,
    BGO **0.88 / 0.77** (= containment²); ≈ 85 % of interacting gammas stay in one
    crystal (CsI), 14 % overspill to a neighbour.

### Navigator — multi-volume transport (Step 3a)

The general walk that carries a photon **across volumes**, switching material at each
boundary: phantom (water, may scatter) → air gap → crystal ring. Built by *reusing* the
single-volume transporter, not duplicating it — three functions, one physics kernel:

- **`sample_interaction`** (`src/sampling.jl`) — the physics at one point (process choice
  + Compton kinematics), record-type agnostic; returns `(process, e_dep, new_dir, new_E)`.
- **`propagate_photon`** (`src/transport.jl`) — the single-volume loop, now returning a
  `Transported{recs, pos, dir, E, escaped}` (the exit state, previously discarded, is what
  lets the navigator chain volumes). Calls the kernel. Single-volume callers read `.recs`.
- **`navigate_photon`** (`src/navigator.jl`) — chains volumes with **no transport loop of
  its own**: `locate` the current volume; in a *leaf* absorber (phantom, ring) transport
  with `propagate_photon` to its surface and carry the exit state on; in the
  non-interacting *air mother* skip straight to `next_boundary`. Output is a tagged stack
  (`NavStep` = `Interaction` + `volume`/`iz`/`iphi`). A backscattered photon crossing the
  bore re-enters the opposite crystal with no special case.
- **`scripts/simulate_phantom.jl`** — back-to-back pairs from the phantom, both
  photons navigated through water → air → ring; CSV tagged with volume, block and a
  per-photon phantom-scatter flag (`--source point|phantom`, `--material`).
- The navigation algorithm (world model, `locate`, `next_boundary`, the leaf/air walk,
  bore re-entry) is documented in `docs/navigation.tex` §4.

**Validated** (`scripts/phantom_effect_on_coincidences.jl`, 20k back-to-back from the
phantom centre, CRYSP1M, CsI — same emission directions air vs phantom):

| per photon | air-only | through the phantom |
|---|---|---|
| unscattered in the phantom | 100 % | **42 %** (Beer–Lambert: exp(−8/10.44) = 46 %) |
| deposited in a crystal | 65 % | 70 % |
| full-energy (≥505 keV, 1 crystal) | 43 % | 18 % |

| per event | air-only | through the phantom |
|---|---|---|
| both photons unscattered | 100 % | **18 %** |
| both full-energy (truth) | 24 % | **4 %** |

**"Full-energy" here is a TRUTH ≥505 keV cut** — at truth level a fully-absorbed photon
deposits exactly 511 keV, so it isolates the full-energy peak, but the cut also rejects
every phantom-scattered photon. So the 4 % is the *unscattered, fully-contained* subset,
**not** the coincidence efficiency: the real selection (Step 4) smears by FWHM(E) and
applies a ±2·FWHM energy window that keeps scattered photons as *scatter* coincidences.
The drop follows `unscat²`: both-full-energy(both-reached) ≈ 0.42² × the air value (8 %
predicted, measured 8 % CsI), because both photons must cross the water unscattered.

The reduction test confirms `navigate_photon` is bit-for-bit `propagate_photon` when a
single absorbing material is crossed (the air leg consumes no randomness).

**Test status:** `Pkg.test` — **202 assertions** pass (foundations, phantom, single
crystal, the `CylShell` shell + the `Sphere` solid, the block/wheel grid, scanner loading,
the air world, the navigator `locate` / `next_boundary` + bore re-entry / reduction /
phantom leg, the emission `Source` + acollinearity, the detector smearing, and the TOML
config). All scripts run.

---

## 3. Next steps

- **Step 2 — detector ring.** *Done* (see above): air world, `CylShell` + the
  (φ, z) block/wheel grid, the `CRYSP_CSI_1M` scanner, and a back-to-back unit test
  reaching the ring through air. Deferred within it: the φ/z plane-crossing
  distances (only needed once dead gaps return).
- **Step 3 — the navigator + coincidences.**
  - *3a — the multi-volume navigator: done* (see above). A photon born in the water
    phantom scatters, exits, crosses the air and reaches the ring, switching material at
    each boundary; validated against Beer–Lambert and the air-only unit test.
  - *3b — coincidences: done.* `scripts/build_coincidences.jl` streams the navigated stack
    into a list-mode coincidence file, tagged **true**/**scatter** from the phantom-scatter
    flag.
- **Step 4 — hits & selection: done.** Hit = first scanner interaction (LOR point) +
  summed crystal energy; detector response (`src/detector.jl`: σ_xyz, FWHM(E)) and the
  energy selection (symmetric window or an `emin_keV` minimum) in `build_coincidences.jl`;
  the two-crystal clean selection (each gamma one block). The whole **analytic-phantom →
  LOR** track (uniform `Source` + acollinearity, `Sphere` solid, TOML-config pipeline,
  per-config + matrix plots, run in parallel) is in `dev/phantom_track_plan.md` (A–F done).
- **Step 5 — randoms.** The time-tag-and-pair pass over the singles (needs real times).
- **Production I/O (Phase G remainder).** The `--singles` reduced stack + HDF5 for 10⁸
  runs (the TOML/naming half is done).
- **Step 6 — detectors.** CsI and BGO done (via the run configs); CsI(Tl), LYSO to add.

Carried-over technical TODOs from the foundations:

- **rotation** in the placement transform, when a volume needs it;
- crystal **XCOM tables** — CsI and BGO done (edges came through numeric); LYSO to
  add, watching the loader's leading-digit filter in case its edge rows are labelled;
- **pair production** in transport, only if runs ever go multi-MeV (correct at 511 keV
  as-is).
