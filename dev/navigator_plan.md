# Plan — Step 3: the multi-volume navigator (+ coincidences)

Self-contained plan so this can be picked up in a fresh session. Status at time of
writing: world (air) + phantom (water) + scanner (CsI `CylShell`) geometry, the
`(φ,z)` block grid, single-volume transport, and the back-to-back unit test (air →
ring, no phantom) are all done and committed; 141 tests pass. The ray–solid math is
in `docs/navigation.tex` (§2–3); the navigation algorithm this plan implements is
documented there too (§4).

## Goal

Transport a photon **across volumes**, switching material at each boundary, in one
history: **phantom (water, scatters) → air gap (straight) → ring (CsI, interacts)**.
This is the piece deferred at every geometry step. It unlocks the phantom-scatter
background and, on top of it, the singles list + same-annihilation coincidences
(true / scatter).

## Architecture — three functions, one physics kernel, no duplication

The design (agreed after review): each concern is written exactly once, and the
multi-volume navigator **reuses** the single-volume transporter rather than
re-implementing its loop.

```
1. sample_interaction   physics at ONE point          (sampling.jl)   ← the kernel
        ▲
2. propagate_photon → Transported{ recs, pos, dir, E, escaped }  (transport.jl)
        ▲                step through ONE solid
        │  (called once per volume; .recs tagged, exit state chained)
3. navigate_photon  → Vector{NavStep}                  (navigator.jl, new)
                         chain solids over the whole Geometry
```

- The **physics** (process choice + Compton kinematics) lives once, in #1.
- The **stepping loop** lives once, in #2.
- The **volume chaining** lives once, in #3 — it has *no* transport loop of its own.

There is **no compatibility wrapper**. The only ripple is that the 7 single-volume
callers of `propagate_photon` read `.recs` off the returned struct (below).

### 1. `sample_interaction` — the physics at a point (`sampling.jl`)

Record-type agnostic. Takes the photon and the already-computed cross sections (so
`sigma_macro` is still called once per step in the loop) and returns the outcome as a
plain tuple — no record, no loop knowledge:

```julia
# returns (process, e_dep, new_dir, new_E);  new_E = 0 if the photon is absorbed
function sample_interaction(E, dir, ΣC, ΣPh, ΣP, rng)
    Σ = ΣC + ΣPh + ΣP
    proc = sample_process(ΣC/Σ, ΣPh/Σ, ΣP/Σ, rng)
    if proc === :compton
        Eprime, cosθ = sample_compton(E, rng)
        ndir = rotate_to_global(scatter_dir(cosθ, rng), dir)
        return (:compton, E - Eprime, ndir, Eprime)   # recoil deposited, photon continues
    else
        return (proc, E, dir, 0.0)                    # photoelectric/pair: full absorption
    end
end
```

This is the current transport.jl:49-56 lifted out, minus record-building and the
below-cut/`break` handling (those are *loop* concerns and differ per record type).

### 2. `propagate_photon` — step through ONE solid (`transport.jl`)

Keeps its name and single-volume loop, calls the kernel, and returns a **rich result**
instead of a bare vector:

```julia
struct Transported
    recs::Vector{Interaction}   # the interaction stack — exactly what callers got before
    pos                         # exit position [cm]
    dir                         # exit direction
    E                           # surviving energy [MeV]
    escaped::Bool               # true = left the volume alive; false = absorbed / below_cut
end

propagate_photon(E, pos, dir, pv, rng; egamma_cut=0.010) -> Transported
```

The loop is unchanged in behaviour (Compton continues, photoelectric/`:below_cut`
stop, an `:escape` `Interaction` is recorded at the exit when the photon leaves alive).
The *only* new thing is that the exit `(pos, dir, E, escaped)` is surfaced rather than
discarded — that is what lets #3 chain volumes.

**Why a struct, not a wrapper.** Returning `Transported` (not `Vector{Interaction}`)
means the 7 existing callers change `propagate_photon(...)` → `propagate_photon(...).recs`.
`.recs` is Julia field access: the same vector as before, pulled out of the struct. No
second function, no compat shim — the exit state is genuinely useful info today's code
drops.

### 3. `navigate_photon` — chain solids over the geometry (`navigator.jl`, new)

No transport loop of its own: it locates the current volume, calls `propagate_photon`
for it, tags the segment, and chains using the exit state.

```julia
function navigate_photon(geom, E0, pos0, dir0, rng; egamma_cut=0.010)::Vector{NavStep}
    pos, dir, E = pos0, dir0, E0
    out = NavStep[]
    nseg = 0
    while true
        vol = locate(geom, pos .+ NAV_EPS .* dir)        # which material now? (nudge past boundary)
        vol === nothing && break                          #   escaped the world
        if vol === geom.world                             # ── air: straight skip, no interaction
            d = next_boundary(geom, vol, pos, dir)
            isfinite(d) || break
            pos = pos .+ d .* dir
        else                                              # ── phantom / ring: a leaf absorber
            seg = propagate_photon(E, pos, dir, vol, rng; egamma_cut)   # reuse #2
            append!(out, tag.(Ref(geom), Ref(vol), seg.recs))          # whole batch = one volume
            seg.escaped || break                          # absorbed / below_cut → done
            pos, dir, E = seg.pos, seg.dir, seg.E         # escaped this volume → carry on
        end
        (nseg += 1) > MAX_SEGMENTS && break               # backstop against a boundary loop
    end
    push!(out, world_escape(pos))                         # final :escape bookkeeping row
    out
end
```

#### Why the air world is special-cased

`propagate_photon` steps to the **current solid's own** `distance_to_exit`, which is
correct only when nothing is embedded inside that solid:

- **phantom, ring** are *leaf* absorbers (nothing nested in them) → `propagate_photon`
  stops exactly at their surface. ✓ Both interact, so transporting them is real work.
- **air world** is the *mother* — it encloses the phantom and ring. `propagate_photon(world)`
  would sail straight to `r = 60` and skip the ring. But air has `Σ = 0` and never
  interacts, so the air leg is not a transport step at all: it is a **straight skip** to
  `next_boundary` (the nearest of phantom-entry / ring-entry / world-exit).

This assumption (one non-interacting mother; phantom/ring are leaves) holds for the
current geometry; a comment flags that a future *nested absorber* would need the
uniform `next_boundary` treatment instead.

#### Bore re-entry falls out for free

A backscattered photon exits the ring's inner wall (a `propagate_photon` escape into
the bore) → lands in air → the air-skip's `next_boundary` returns the **far** wall →
`propagate_photon` runs again on the opposite crystal. The overspill-across-bore the
single-volume unit test cannot follow is handled with no special code.

### `locate(geom, p) -> PhysicalVolume | nothing`

Innermost daughter first, mother last:
```julia
is_inside(geom.phantom, p)               && return geom.phantom
geom.scanner !== nothing &&
    is_inside(geom.scanner.volume, p)    && return geom.scanner.volume
is_inside(geom.world, p)                 && return geom.world
return nothing   # escaped the world
```

### `next_boundary(geom, current, pos, dir) -> Float64`

Nearest distance at which volume membership changes (exit `current`, or enter another):
```julia
d = distance_to_exit(pos, dir, current)
for v in others(geom, current)            # the volumes ≠ current (skip scanner if nothing)
    d = min(d, distance_to_entry(pos, dir, v))
end
return d
```
Correct because `distance_to_entry` returns `Inf` for a point already inside, so being
inside the mother never produces a spurious candidate. Used by the air-skip in #3 (and
available for the bore-re-entry test).

### Record: `NavStep` ( = `Interaction` + location tags )

```julia
struct NavStep
    hit::Interaction      # the 6 physics fields, unchanged (composition, not duplication)
    volume::Symbol        # :phantom | :air | :scanner
    iz::Int               # block, or -1 off the scanner
    iphi::Int
end
```
`tag(geom, vol, interaction)` builds it: stamps `volume` from `vol` and fills `(iz,iphi)`
via `block_index(geom.scanner, pos)` only when `vol === geom.scanner.volume`, else
`(-1,-1)`. Single source of truth for the physics fields — if `Interaction` gains a
field (e.g. time) `NavStep` inherits it free.

### Constants

`NAV_EPS ≈ 1e-7 cm` (1 nm) — the nudge to re-locate past a just-crossed boundary; the
distance functions' `SURFACE_EPS = 1e-10` already rejects the surface itself.
`MAX_SEGMENTS` — a stuck-on-boundary backstop.

## The escape rows (bookkeeping)

Escape rows fall out naturally and form a faithful per-volume trace:
- `propagate_photon(phantom)` ends with `:escape` at the phantom surface (tagged
  `:phantom`) — the photon left the patient.
- `propagate_photon(ring)` ends with `:escape` (overspill out / backscatter) **or** a
  terminal `:photoelectric` / `:below_cut` (absorbed, no escape row), correctly.
- a final `:escape` row when the photon leaves the world (`world_escape` above).

## Source from the phantom

Annihilations now originate **in the phantom**, not a bare point. For Step 3:
- start with the existing back-to-back point source **at the origin** — but now it flies
  through the water phantom first (validation: phantom scatter appears);
- then a phantom-distributed source: a point uniform in the phantom cylinder +
  back-to-back isotropic pair;
- the scenario source (`emitters.csv`) is later (pipeline step 1).

## Scope (split)

- **3a (do first):** the navigator library (`sample_interaction`, the `Transported`
  rework of `propagate_photon`, `locate`, `next_boundary`, `navigate_photon`, `NavStep`,
  `tag`), the 7-caller `.recs` rework, reworked tests, and a validation script
  `scripts/navigate_back_to_back.jl` (emit back-to-back from the origin / phantom,
  navigate both photons, write the stack tagged with volume + block + a phantom-scatter
  flag).
- **3b (next):** phantom-distributed source + same-annihilation coincidence output — a
  coincidence = both photons produce a scanner deposit; tag **true** (neither scattered
  in the phantom) / **scatter** (≥1 did). The detailed hit formation (first-interaction
  point, smear, energy window, two-block selection) is **Step 4**, not here.

## Files to touch

- `src/sampling.jl` — add `sample_interaction` (the shared kernel).
- `src/transport.jl` — `propagate_photon` returns `Transported`, calls the kernel.
- `src/navigator.jl` (new) — `locate`, `next_boundary`, `navigate_photon`, `NavStep`,
  `tag`, `NAV_EPS`, `MAX_SEGMENTS`.
- `src/PTCryspMC.jl` — `include("navigator.jl")` + exports (`navigate_photon`,
  `NavStep`, `locate`, `next_boundary`, `Transported`).
- 4 scripts (`shoot_back_to_back_511_keV_gammas.jl`, `shoot_gammas_to_crystal.jl`,
  `propagate_gammas_in_phantom.jl`, `bench_back_to_back.jl`) — `propagate_photon(...)`
  → `propagate_photon(...).recs`.
- `test/runtests.jl` — the two air tests (287, 295) → `.recs`; **replace** the
  source→ring test (342, now a subset of `navigate_photon`) with a navigator testset.
- `scripts/navigate_back_to_back.jl` (new) — validation/driver.
- Optional later: a plotter; update `dev_steps.md` + `pet_simulation.tex`.

## Validation / tests (reworked, no compat shims)

- **Reduction**: a photon fired through air straight at the ring (missing the phantom)
  gives the **same** scanner deposits as a single-volume `propagate_photon(...).recs`
  into the scanner — `navigate_photon` reduces to single-volume when one material is
  crossed. (Replaces the old source→ring test.)
- **Phantom leg**: a photon from the phantom centre records **water (Compton)**
  interactions before reaching the ring; the phantom-scatter flag is set iff it
  Compton'd in water.
- **Energy conservation** across volumes (Σ e_dep ≤ E0; = E0 if fully absorbed).
- **Bore re-entry capability**: from the bore, `next_boundary` returns the **far** wall
  (backscatter can reach the opposite crystal).
- **Phantom degrades coincidences**: back-to-back through water gives a **lower**
  clean-coincidence fraction than the air-only test (0.29 CsI) — the phantom-scatter
  loss appears.

## Risks / watch-outs

- **Boundary epsilon**: too small a `NAV_EPS` → re-locate lands on the boundary
  (ambiguous, inclusive `is_inside` true for both); too large → a visible position error.
  ~1e-7 cm (1 nm) is safe.
- **Segment loop guard**: if the air-skip / leaf hand-off ever stalls on a boundary,
  `MAX_SEGMENTS` backstops it.
- **Air `Σ=0`** → handled by the dedicated skip (not by `propagate_photon`).
- Don't change `propagate_photon`'s *behaviour* — only its return type; the reworked
  tests must still see identical physics.

## Quick reference (numbers)

CRYSP1M ring: Ri=38.7, wall=3.7, Ro=42.4, H=51.2 cm; 48 φ × 20 z = 960. Phantom: water
cylinder r=8, H=8 cm. World: air cylinder r=60, H=60 cm. Air-only unit test: acceptance
79.5%, clean-coincidence CsI 0.29 / BGO 0.77.
