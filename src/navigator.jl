# Multi-volume navigation: carry a photon across volumes (water phantom → air gap →
# CsI ring), switching the material at each boundary, by reusing the single-volume
# `propagate_photon` rather than re-implementing its loop. See docs/navigation.tex §4
# and dev/navigator_plan.md.
#
# World model: one non-interacting Air mother enclosing two radially disjoint leaf
# absorbers (phantom, ring). A photon inside a leaf is transported by `propagate_photon`,
# which stops at that leaf's own surface; the air mother — being Σ=0 — is never
# transported, only skipped to the next interface (`next_boundary`).

const NAV_EPS      = 1e-7    # [cm] nudge to re-locate past a just-crossed boundary (~1 nm)
const MAX_SEGMENTS = 1000    # backstop against a photon stuck on a boundary

"""
    NavStep

One recorded interaction of a navigated photon: the physics (`hit::Interaction`) plus
where it happened — `volume` (`:phantom` | `:air` | `:scanner`) and, in the scanner,
the readout block `(iz, iphi)` (`-1` off the scanner). Composition keeps a single
source of truth for the physics fields (access them as `step.hit.e_dep`, …).
"""
struct NavStep
    hit::Interaction
    volume::Symbol
    iz::Int
    iphi::Int
end

"""
    locate(geom, p) -> PhysicalVolume | Nothing

The volume whose material applies at world-frame point `p`: innermost daughter first
(phantom, then the scanner ring), the Air mother last, or `nothing` if `p` has escaped
the world.
"""
function locate(geom::Geometry, p)
    is_inside(geom.phantom, p) && return geom.phantom
    geom.scanner !== nothing && is_inside(geom.scanner.volume, p) && return geom.scanner.volume
    is_inside(geom.world, p) && return geom.world
    nothing
end

"""
    next_boundary(geom, current, pos, dir) -> Float64

Nearest distance at which volume membership changes: leave `current`, or enter another
volume. `distance_to_entry` is `Inf` for a volume already containing `pos`, so the
enclosing mother never yields a spurious candidate. Used for the air traverse — the
mother is non-interacting, so its step is a straight skip to here (the nearest of
phantom-entry, ring-entry and world-exit).
"""
function next_boundary(geom::Geometry, current, pos, dir)::Float64
    d = distance_to_exit(pos, dir, current)
    current === geom.phantom || (d = min(d, distance_to_entry(pos, dir, geom.phantom)))
    current === geom.world   || (d = min(d, distance_to_entry(pos, dir, geom.world)))
    if geom.scanner !== nothing && current !== geom.scanner.volume
        d = min(d, distance_to_entry(pos, dir, geom.scanner.volume))
    end
    d
end

# Tag a single-volume interaction with its location: the block (iz, iphi) for a scanner
# deposit, (-1, -1) otherwise.
@inline function _tag(geom::Geometry, tagvol::Symbol, r::Interaction)::NavStep
    if tagvol === :scanner
        iφ, iz = block_index(geom.scanner, (r.x, r.y, r.z))
        NavStep(r, :scanner, iz, iφ)
    else
        NavStep(r, tagvol, -1, -1)
    end
end

"""
    navigate_photon(geom, E0, pos0, dir0, rng; egamma_cut=0.010) -> Vector{NavStep}

Transport one photon of energy `E0` [MeV] from `pos0` along `dir0` across the whole
geometry, switching material at every boundary. The walk reuses `propagate_photon` for
each leaf volume (phantom, ring) and skips straight through the Air mother to the next
interface; it has no transport loop of its own. Returns the tagged interaction stack
(each `NavStep` carries its volume and, in the ring, its block). A backscattered photon
crossing the bore re-enters the opposite crystal with no special case.
"""
function navigate_photon(geom::Geometry, E0::Real, pos0, dir0, rng::AbstractRNG;
                         egamma_cut::Float64=0.010)::Vector{NavStep}
    pos = (Float64(pos0[1]), Float64(pos0[2]), Float64(pos0[3]))
    dir = (Float64(dir0[1]), Float64(dir0[2]), Float64(dir0[3]))
    n   = sqrt(dir[1]^2 + dir[2]^2 + dir[3]^2)
    dir = (dir[1]/n, dir[2]/n, dir[3]/n)
    E   = Float64(E0)

    out = NavStep[]
    escaped_world = true
    nseg = 0
    while true
        probe = (pos[1] + NAV_EPS*dir[1], pos[2] + NAV_EPS*dir[2], pos[3] + NAV_EPS*dir[3])
        vol = locate(geom, probe)
        vol === nothing && break                       # left the world

        if vol === geom.world                          # ── air: straight skip, no interaction
            d = next_boundary(geom, vol, pos, dir)
            isfinite(d) || (escaped_world = false; break)
            pos = (pos[1] + d*dir[1], pos[2] + d*dir[2], pos[3] + d*dir[3])
        else                                           # ── phantom / ring: a leaf absorber
            tagvol = vol === geom.phantom ? :phantom : :scanner
            seg = propagate_photon(E, pos, dir, vol, rng; egamma_cut=egamma_cut)
            for r in seg.recs
                push!(out, _tag(geom, tagvol, r))
            end
            if !seg.escaped                            # absorbed / below cut → history ends
                escaped_world = false
                break
            end
            pos, dir, E = seg.pos, seg.dir, seg.E      # escaped this volume → carry on
        end

        nseg += 1
        nseg > MAX_SEGMENTS && (escaped_world = false; break)
    end

    # Bookkeeping: a final :escape row only when the photon actually leaves the world
    # (not when it was absorbed in a leaf, where the terminal record is the absorption).
    escaped_world && push!(out, NavStep(Interaction(pos[1], pos[2], pos[3], E, 0.0, :escape),
                                        :air, -1, -1))
    out
end
