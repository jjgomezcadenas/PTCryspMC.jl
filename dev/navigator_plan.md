# Plan — Step 3: the multi-volume navigator (+ coincidences)

Self-contained plan so this can be picked up in a fresh session. Status at time of
writing: world (air) + phantom (water) + scanner (CsI `CylShell`) geometry, the
`(φ,z)` block grid, single-volume transport, and the back-to-back unit test (air →
ring, no phantom) are all done and committed; 141 tests pass. The math for ray–solid
distances is in `docs/navigation.tex`.

## Goal

Transport a photon **across volumes**, switching material at each boundary, in one
history: **phantom (water, scatters) → air gap (straight) → ring (CsI, interacts)**.
This is the piece deferred at every geometry step. It unlocks the phantom-scatter
background and, on top of it, the singles list + same-annihilation coincidences
(true / scatter).

## What already exists (reuse, don't rebuild)

- `Geometry{world::PhysicalVolume(Air), phantom::PhysicalVolume(Water), scanner::Union{Scanner,Nothing}}`
  in `src/geometry.jl`.
- Per-solid `is_inside`, `distance_to_entry`, `distance_to_exit` (Cylinder, Box,
  CylShell) and their `PhysicalVolume` delegations (world→local via `_to_local`).
- `propagate_photon(E0, pos, dir, pv, rng; egamma_cut)` — single-volume loop,
  returns `Vector{Interaction}` (x,y,z,e_in,e_dep,process). **Keep it** (used by the
  crystal/box studies and the air-only unit test).
- Samplers in `sampling.jl` (`sample_distance`, `sample_process`, `sample_compton`,
  `rotate_to_global`); `sigma_macro`; `block_index(scanner, p)`/`block_id`.
- `egamma_cut` default 0.010 MeV (10 keV, XCOM min). `SURFACE_EPS=1e-10`,
  `PARALLEL_EPS=1e-20` in geometry.jl.

## The world model

Air world (mother) [0,60] cm contains the water phantom [0,8] and the CsI shell
[38.7,42.4] cm. Daughters do **not** overlap each other (phantom and ring are
radially separate). Membership of a point: phantom if inside phantom; else scanner
if inside scanner; else world (air) if inside world; else escaped.

## Design

New file `src/navigator.jl` (include + export from `PTCryspMC.jl`). Three pieces:

### 1. `locate(geom, p) -> PhysicalVolume | nothing`
Innermost daughter first, mother last:
```julia
is_inside(geom.phantom, p)               && return geom.phantom
geom.scanner !== nothing &&
    is_inside(geom.scanner.volume, p)    && return geom.scanner.volume
is_inside(geom.world, p)                 && return geom.world
return nothing   # escaped the world
```
Returns the volume whose **material** to use. (Block tagging uses `geom.scanner`
separately, only for scanner deposits.)

### 2. `next_boundary(geom, current, pos, dir) -> Float64`
The nearest distance at which volume membership changes: leave the current volume,
or enter another one.
```julia
d = distance_to_exit(pos, dir, current)
for v in (geom.phantom, geom.scanner.volume, geom.world)   # skip current; skip scanner if nothing
    v === current && continue
    d = min(d, distance_to_entry(pos, dir, v))
end
return d
```
Notes / why it's correct:
- From inside the phantom (⊂ world): `distance_to_entry(world)=Inf` (already inside),
  so the candidate is `exit(phantom)` vs `entry(scanner)` → leaves into air first. ✓
- From the air world heading at the ring: `entry(scanner)` (near face) vs
  `entry(phantom)` vs `exit(world)` (far). ✓
- From the CsI shell, a backscattered photon exits the **inner** wall into the bore
  (`exit(scanner)` = inner wall) → air; then from the bore `entry(scanner)` = the
  **far** wall, so it can **re-enter the opposite crystal**. This is the
  overspill-across-bore the single-volume unit test could not follow — the navigator
  gets it for free. (Tested as a capability.)

### 3. `navigate_photon(geom, E0, pos0, dir0, rng; egamma_cut) -> Vector{NavStep}`
Generalises `propagate_photon`'s loop with material switching:
```
pos, dir, E = ...
while true
    vol = locate(geom, pos + NAV_EPS*dir)   # nudge to disambiguate the boundary
    vol === nothing && break                  # escaped
    mat = material(vol)
    Σ = sum(sigma_macro(mat, E)); s = Σ>0 ? sample_distance(Σ,rng) : Inf
    d_b = next_boundary(geom, vol, pos, dir)
    if s >= d_b                               # cross a boundary, no interaction
        pos += d_b*dir
        continue                              # re-locate next iteration (nudged)
    else                                      # interact
        pos += s*dir
        Compton (deposit recoil, new dir, E') / photoelectric (deposit all, stop) /
        below_cut (deposit E, stop)
        record NavStep(pos, E, e_dep, process, volume_tag, iz, iphi)
        E < egamma_cut && break
    end
end
```
- `NAV_EPS` ≈ `1e-7` cm (boundary-crossing nudge for robust re-location; the
  `SURFACE_EPS` in the distance fns already rejects the just-crossed surface via
  `t>SURFACE_EPS`). Pick/validate during impl.
- The air legs are a single straight step (Σ=0 → s=Inf → cross to the next boundary).

### Record: `NavStep`
The navigator needs, per interaction, the **volume** (so we can tell a phantom
scatter from a detector deposit) and, for scanner deposits, the **block**. Decide at
impl time between:
- (a) a new `struct NavStep` (x,y,z,e_in,e_dep,process,volume::Symbol,iz::Int,iphi::Int;
  iz/iphi = -1 off the scanner), or
- (b) extend `Interaction` with the extra fields (touches single-volume transport).
Lean **(a)** to keep `Interaction`/`propagate_photon` untouched. `volume` ∈
{:phantom,:air,:scanner}. Compute `(iz,iphi)` via `block_index(geom.scanner, pos)`
only when `vol === geom.scanner.volume`.

## Source from the phantom

Annihilations now originate **in the phantom**, not a bare point. For Step 3:
- start with the existing back-to-back point source **at the origin** — but now it
  flies through the water phantom first (validation: phantom scatter appears).
- then a phantom-distributed source: a point uniform in the phantom cylinder +
  back-to-back isotropic pair.
- the scenario source (`emitters.csv`) is later (pipeline step 1).
Move `emit_pair` toward a small `Source` (we flagged this); keep it script-local /
minimal for now, sharing the `(pos, dir1, dir2, E)` shape.

## Scope (split)

- **3a (do first):** the navigator library (`locate`, `next_boundary`,
  `navigate_photon`, `NavStep`), exports, tests, and a validation script
  `scripts/navigate_back_to_back.jl` (emit back-to-back from the origin / phantom,
  navigate both photons, write the stack tagged with volume + block + a
  phantom-scatter flag).
- **3b (next):** phantom-distributed source + same-annihilation coincidence output —
  a coincidence = both photons produce a scanner deposit; tag **true** (neither
  scattered in the phantom) / **scatter** (≥1 did). The detailed hit formation
  (first-interaction point, smear, energy window, the two-opposite-block clean
  selection) is **Step 4**, not here.

## Files to touch

- `src/navigator.jl` (new) — `locate`, `next_boundary`, `navigate_photon`, `NavStep`.
- `src/PTCryspMC.jl` — `include("navigator.jl")` + exports.
- `test/runtests.jl` — new testset (below).
- `scripts/navigate_back_to_back.jl` (new) — validation/driver.
- Optional later: a plotter; update `dev_steps.md` + `pet_simulation.tex`.

## Validation / tests

- A photon launched in air aimed straight at the ring (missing the phantom) gives
  the **same** result as the single-volume `propagate_photon` into the scanner
  (navigator reduces to single-volume when only one material is crossed).
- A photon from the phantom centre records **water (Compton) interactions** before
  reaching the ring; the phantom-scatter flag is set iff it Compton'd in water.
- **Energy conservation** across volumes (Σ e_dep ≤ E0; = E0 if fully absorbed).
- The **back-to-back through water** clean-coincidence fraction is **lower** than the
  air-only unit test (0.29 CsI) — the phantom-scatter loss/background appears.
- Capability: a backscattered photon crossing the bore **re-enters** the opposite
  crystal (construct a geometry-level case and check `next_boundary` returns the far
  wall from the bore).

## Risks / watch-outs

- **Boundary epsilon**: too small a `NAV_EPS` → re-locate lands on the boundary
  (ambiguous, inclusive `is_inside` true for both); too large → a visible position
  error. ~1e-7 cm (1 nm) is safe.
- **Infinite loop guard**: if `next_boundary` ever returns ~0 repeatedly (stuck on a
  boundary), add a max-steps cap per photon as a backstop.
- **Air `Σ=0`** → `s=Inf`; the loop must cross to the boundary (works: `s>=d_b`).
- Don't break `propagate_photon` / the existing 141 tests.

## Quick reference (numbers)

CRYSP1M ring: Ri=38.7, wall=3.7, Ro=42.4, H=51.2 cm; 48 φ × 20 z = 960. Phantom:
water cylinder r=8, H=8 cm. World: air cylinder r=60, H=60 cm. Air-only unit test:
acceptance 79.5%, clean-coincidence CsI 0.29 / BGO 0.77.
