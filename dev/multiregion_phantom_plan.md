# Multi-region (non-uniform) phantom — implementation plan

**Status: deferred, fully scoped. Not yet built.** This is the "MIRD multi-region" extension
long noted in `dev/status.md`. The single-region phantom path is built and validated
(`load_phantom_regions` in `src/scenario.jl`, `Ellipsoid` in `src/geometry.jl`). This document is
a self-contained brief for a future instance to implement multi-region support **without
re-investigating** — the data, the design, the prerequisites, and the tests are all here.

Read `dev/api_plan.md` first for the surrounding API-source work; this is one deferred piece of it.

---

## Why (motivation)

The uniform head (`uniform_headep_sobp_1e8`) fills the whole head envelope with one material
(brain). The realistic head phantoms layer **soft tissue / brain / cortical bone**, and bone
attenuates 511 keV photons **~1.7× more than brain** (μ 0.171 vs 0.099 cm⁻¹). That skull-shell
contrast measurably shapes both the transported 511 keV flux and the attenuation-correction map a
reconstruction needs, so a faithful range-verification study on a head wants the layered phantom,
not the uniform approximation.

## The data (from ptcryspg4, verified)

Two multi-region scenarios, both **nested ellipsoids resolved by priority** (lowest `priority`
number wins where regions overlap — the documented "brain-first carves the skull shell" model),
all `ellipsoid`, all Euler = 0. Confirmed geometrically nested in priority order, with the
outermost region (scalp) bounding all others.

`.../ptcryspg4/data/runs/headep_sobp_1e7/phantom_regions.csv` — 4 regions, 3 materials:

| priority | region | material | semi-axes a,b,c mm | centre mm |
|---|---|---|---|---|
| 0 | tumour | G4_TISSUE_SOFT_ICRP | 18,18,20 | 0,0,-25 |
| 1 | brain  | G4_BRAIN_ICRP        | 60,65,90 | 0,-40,0 |
| 2 | skull  | G4_BONE_CORTICAL_ICRP| 68,83,98 | 0,-30,0 |
| 3 | scalp  | G4_TISSUE_SOFT_ICRP  | 72,87,102| 0,-30,0 |

`.../ptcryspg4/data/runs/mird_head_sobp_1e7/phantom_regions.csv` — 3 regions, 3 materials:

| priority | region | material | semi-axes a,b,c mm | centre mm |
|---|---|---|---|---|
| 0 | brain | G4_BRAIN_ICRP         | 65,90,60 | 0,0,0 |
| 1 | skull | G4_BONE_CORTICAL_ICRP | 83,98,68 | -10,0,0 |
| 2 | scalp | G4_TISSUE_SOFT_ICRP   | 87,102,72| -10,0,0 |

Material lookup at a point: the first (lowest-priority) region containing it, else air —
`((x-cx)/a)² + ((y-cy)/b)² + ((z-cz)/c)² ≤ 1` per ellipsoid.

---

## Prerequisites

### 1. Two new materials + their XCOM tables (BLOCKER)

`G4_BRAIN_ICRP` is done (`data/materials.json`, `data/xcom_brain.csv`). Needed:

| material | ρ g/cm³ | μ(511) cm⁻¹ (validation target) | composition source |
|---|---|---|---|
| G4_BONE_CORTICAL_ICRP | 1.92 | 0.170555 (μ/ρ 0.088831) | frozen scenario dir |
| G4_TISSUE_SOFT_ICRP   | 1.03 | 0.098010 (μ/ρ 0.095155) | frozen scenario dir |

Compositions (mass fractions), from `ptcrysp-scenarios/.../phantom_material_g4_*_icrp.csv`:
- **bone_cortical**: H .034, C .155, N .042, O .435, Na .001, Mg .002, P .103, S .003, Ca .225
- **tissue_soft**:   H .105, C .256, N .027, O .602, Na .001, P .002, S .003, Cl .002, K .002

**XCOM tables** `data/xcom_bone_cortical.csv`, `data/xcom_tissue_soft.csv` are NOT present. Fetch
from NIST XCOM as the brain/water tables were (the project convention), or generate the mixture
from the compositions. Note: Ca (K-edge 4.0 keV), P (2.1), K (3.6), Cl (2.8) all sit **below** the
10 keV grid floor, so there is **no in-range K-edge** — the current `load_xcom` parser handles
these fine (does not trigger the deferred K-edge hardening, `dev/status.md` #2).

`materials.json` entries mirror the brain entry (density, components, xcom path; no scintillator
props). Add a μ-validation test per material against the targets above (as done for brain in the
"materials loading" testset), tolerance ~1%.

### 2. The single-region foundation (done)

`Ellipsoid` solid, `load_phantom_regions` (single-region, errors on >1), the placement/units
convention (mm→cm, centre → `PhysicalVolume.position`). Multi-region builds directly on these.

---

## Design (additive — do NOT disturb the single-region hot path)

The navigator today has **one phantom leaf, one material** (`_leaf_reduce(geom.phantom, …)`
transports to exit using `material(geom.phantom)`). Multi-region means the material **changes
mid-flight** across region boundaries. Keep the existing single-`PhysicalVolume` path exactly as
is (all current configs + tests rely on it); add a parallel multi-region path selected by dispatch.

### `MultiRegionPhantom` (new type, `src/geometry.jl`)

```
struct MultiRegionPhantom
    regions::Vector{PhysicalVolume{Ellipsoid}}   # PRIORITY ORDER (index 1 = priority 0, wins)
    bound::Int                                    # index of the bounding region (contains all others)
end
```
All regions are `Ellipsoid`, so the vector is **concretely typed → no boxing**, preserving the
alloc-free guarantee. **Design constraint:** require a bounding region (the outermost, largest —
scalp here) that geometrically contains every other region; assert this at load. The bounding
region then behaves exactly like today's single phantom solid for the top-level navigator
(`is_inside`, `distance_to_entry`), confining all new logic to the *internal* stepping.

Accessors:
- `material_at(ph, p)` → the material of the first region (priority order) with `is_inside(region,p)`;
  `nothing` if none (outside the union → treat as exited).
- `is_inside(ph::MultiRegionPhantom, p)` → `is_inside(ph.regions[ph.bound], p)` (the bounding region).
- `distance_to_entry/exit(pos,dir, ph)` → delegate to the bounding region (top-level use only).

### Region-aware transport step (`src/navigator.jl`)

A `_leaf_reduce` variant for the phantom when it is a `MultiRegionPhantom` (INSCAN=false always —
the phantom never has blocks). Per step, inside the bounding region:
1. `mat = material_at(ph, pos + NAV_EPS·dir)`; if `nothing` → left the union → `escaped=true`, break.
2. Sample free path `s` with `mat`.
3. `d_bound` = nearest crossing of ANY region surface ahead (the material-winner can only change at
   a region surface): `min` over regions of `distance_to_exit` (if it contains pos) / `distance_to_entry`
   (if not), positive only.
4. If `s < d_bound`: interact at `pos + s·dir` — same Compton/photoelectric/pair/below-cut logic as
   the current phantom branch of `_leaf_reduce` (increment `nscat` on a real interaction, respecting
   the `:below_cut` rule from the reducer fix). Else advance to `pos + (d_bound+NAV_EPS)·dir` and
   re-loop (material re-evaluated at the boundary).

This mirrors the top-level world→phantom→ring navigation one level down. For 3–4 regions the
per-step region scan is a handful of `is_inside`/`distance` calls — cheap.

### Navigator dispatch

`navigate_single_photons` and `locate_tag`/`next_boundary` dispatch on `typeof(geom.phantom)`:
a `PhysicalVolume` → today's path unchanged; a `MultiRegionPhantom` → the variant above.
`Geometry{W,P}` is already parametric on the phantom type, so both stay concretely typed.
**Extend the dev shell `navigate_photon` too** (same region-aware stepping, recording every
NavStep) so the two-shell bit-equality contract (`_reduce_steps(navigate_photon) ==
navigate_single_photons`) still holds for multi-region — do not let the shells diverge.

### Loader (`src/scenario.jl`)

Lift the current "error on >1 region" in `load_phantom_regions`: build a `MultiRegionPhantom` from
all rows (sorted by `priority`), each row → `PhysicalVolume{Ellipsoid}` as today; assert Euler = 0
per row; assert a bounding region exists (all others ⊂ it) and record its index. Single-region
returns the existing single `PhysicalVolume` (back-compat), N>1 returns a `MultiRegionPhantom`.

---

## Tests & validation

- **Priority lookup**: `material_at` on the real headep/mird stacks returns tumour/brain/skull/scalp
  materials at representative points; air outside the bounding region.
- **Attenuation through bone**: a 511 keV pencil beam along a chord that crosses the skull shell
  shows the extra attenuation vs the same chord in a uniform-brain head (Beer–Lambert with the
  per-region μ over the path segments) — the physics reason for the whole feature.
- **Zero allocation**: the multi-region `navigate_single_photons` allocates nothing per photon
  (the `@allocated` slope test, as for the single-region path).
- **Two-shell equality**: `_reduce_steps(navigate_photon(...)) == navigate_single_photons(...)` on
  the layered head over many seeds (the drift guard, extended to multi-region).
- **Per-material μ**: `sigma_macro` reproduces the bone/tissue meta μ to ~1% (prerequisite tests).
- **Loader**: real headep (4) and mird (3) files build a `MultiRegionPhantom` with the right
  region count, priority order, and bounding index; a non-bounding-region file (no region contains
  the rest) errors.

---

## Still deferred after this (out of scope even for multi-region)

- **Region rotation** (nonzero Euler): still asserted zero — the rotation transform in
  `PhysicalVolume._to_local` is a separate deferred item (`dev/status.md`).
- **Non-ellipsoid multi-region** (mixed solid types in one phantom): the homogeneous
  `Vector{PhysicalVolume{Ellipsoid}}` assumption would box; revisit only if a scenario needs it.
- **Overlapping non-nested regions with no single bounding volume**: the bounding-region constraint
  rules this out; the priority stepping itself is general, but the top-level integration assumes a
  bound. Revisit if a scenario violates it.

---

## Instructions for a future instance

1. Read `dev/api_plan.md` (surrounding work, now built) and this file. The single-region foundation
   (`Ellipsoid`, `load_phantom_regions`, `G4_BRAIN_ICRP`) is committed on `main`.
2. **Confirm the prerequisite** with the user: are `data/xcom_bone_cortical.csv` and
   `data/xcom_tissue_soft.csv` present (NIST-fetched)? If not, either they fetch them or generate
   from the compositions above. Do not proceed without the tables.
3. Add the two `materials.json` entries + μ-validation tests (targets above). Land + test first.
4. Build `MultiRegionPhantom` + `material_at` + accessors; then the region-aware `_leaf_reduce`
   variant in **both** navigator shells; then the loader lift. Keep the single-region path and its
   tests untouched — verify the full existing suite stays green at each step.
5. Add the tests above. Validate the bone-attenuation physics explicitly (it is the point).
6. Add `runs/*.toml` + driver support for a multi-region API scenario once the source path
   (`dev/api_plan.md` steps 4–7) exists.
7. Update `dev/status.md` (move multi-region from deferred to done) and `dev/api_plan.md` (fold in).
