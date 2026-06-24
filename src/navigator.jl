# Multi-volume navigation: carry a photon across volumes (water phantom → air gap →
# CsI ring), switching the material at each boundary, by reusing the single-volume
# `propagate_photon` rather than re-implementing its loop. See docs/PTCryspMC_phys.tex §4
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

("navigate_photon_full_stack".) Records EVERY interaction → a `Vector{NavStep}`. Allocates
(the vector here, plus `propagate_photon`'s `recs` per leaf), so it is the dev / detailed-plot
/ CNN-truth path used by `simulate_phantom.jl`. For the production singles path use the
allocation-free `navigate_single_photons` below.

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

# The singles summary one photon reduces to: did it reach the ring, the first crystal
# interaction (the LOR point + its block), the summed crystal energy, the distinct-block
# count (1 = contained, >1 = overspill), and whether it interacted in the phantom. Positions
# in cm (as in `Interaction`); the driver converts to mm on write.
const _EMPTY_SINGLE = (reached=false, x=0.0, y=0.0, z=0.0, iz=-1, iphi=-1,
                       e=0.0, nblocks=0, phscat=false)

"""
    _reduce_steps(steps) -> NamedTuple

Reduce a full-stack `Vector{NavStep}` (the `navigate_photon` output) to the singles summary,
by the same rule `build_coincidences` uses: the first scanner deposit (process ≠ `:escape`)
fixes the LOR point and block; later scanner deposits add energy and bump the distinct-block
count when the block changes; any phantom interaction sets `phscat`. The single source of
truth for the reduction — `navigate_single_photons` must return exactly this on the same
history (asserted in the tests), and the future singles-reader reuses it.
"""
function _reduce_steps(steps)
    reached = false
    hx = 0.0; hy = 0.0; hz = 0.0; hiz = -1; hiphi = -1
    esum = 0.0; nblk = 0; pbz = -2; pbphi = -2
    phscat = false
    for st in steps
        st.hit.process === :escape && continue
        if st.volume === :scanner
            if !reached
                reached = true; hx = st.hit.x; hy = st.hit.y; hz = st.hit.z
                hiz = st.iz; hiphi = st.iphi; nblk = 1; pbz = st.iz; pbphi = st.iphi
            elseif st.iz != pbz || st.iphi != pbphi
                nblk += 1; pbz = st.iz; pbphi = st.iphi
            end
            esum += st.hit.e_dep
        elseif st.volume === :phantom
            phscat = true
        end
    end
    (reached=reached, x=hx, y=hy, z=hz, iz=hiz, iphi=hiphi, e=esum, nblocks=nblk, phscat=phscat)
end

# Like `locate` but returns a Symbol TAG (:phantom / :scanner / :world / :none) instead of
# the volume object. `locate` returns a `PhysicalVolume`, which is a non-isbits immutable
# (it points at a String + Material), so returning it through a Union return BOXES it (~144 B
# per call); at ~3 segments/photon that is the ~455 B/photon the allocation profiler traced
# to here. The tag is a Symbol (interned → no allocation); the caller then transports through
# the concrete `geom` field directly, never holding the volume as a Union.
@inline function locate_tag(geom::Geometry, p)::Symbol
    is_inside(geom.phantom, p) && return :phantom
    geom.scanner !== nothing && is_inside(geom.scanner.volume, p) && return :scanner
    is_inside(geom.world, p) && return :world
    :none
end

# Transport one photon through a single concretely-typed leaf `pv`, folding each interaction
# into the running singles accumulators (mirrors `propagate_photon`'s step loop). `INSCAN`
# (a `Val`) is compile-time, so the scanner-only block/overspill bookkeeping is dead-code
# eliminated for the phantom call (where `sc` is `nothing`). Called with the concrete
# `geom.phantom` / `geom.scanner.volume`, so every `material`/`distance_to_exit`/`block_index`
# is statically dispatched — no Union, no boxing. Returns the updated transport + accumulator
# state as one isbits tuple.
@inline function _leaf_reduce(pv, ::Val{INSCAN}, sc, pos, dir, E::Float64, egamma_cut::Float64,
                              reached::Bool, hx::Float64, hy::Float64, hz::Float64,
                              hiz::Int, hiphi::Int, esum::Float64, nblk::Int,
                              pbz::Int, pbphi::Int, phscat::Bool, rng::AbstractRNG) where {INSCAN}
    nn  = sqrt(dir[1]^2 + dir[2]^2 + dir[3]^2)          # mirror propagate_photon's entry norm
    dir = (dir[1]/nn, dir[2]/nn, dir[3]/nn)
    mat = material(pv)
    escaped = false
    while true
        ΣC, ΣPh, ΣP = sigma_macro(mat, E)
        Σ = ΣC + ΣPh + ΣP
        s = Σ > 0.0 ? sample_distance(Σ, rng) : Inf
        d_exit = distance_to_exit(pos, dir, pv)
        if s >= d_exit                                 # leaves this volume alive
            pos = (pos[1] + d_exit*dir[1], pos[2] + d_exit*dir[2], pos[3] + d_exit*dir[3])
            escaped = true
            break
        end
        pos = (pos[1] + s*dir[1], pos[2] + s*dir[2], pos[3] + s*dir[3])
        is_comp, e_dep, ndir, nE = sample_interaction_t(E, dir, ΣC, ΣPh, ΣP, rng)
        if INSCAN
            iφ, iz = block_index(sc, pos)
            if !reached
                reached = true; hx = pos[1]; hy = pos[2]; hz = pos[3]
                hiz = iz; hiphi = iφ; nblk = 1; pbz = iz; pbphi = iφ
            elseif iz != pbz || iφ != pbphi
                nblk += 1; pbz = iz; pbphi = iφ
            end
            esum += e_dep
        else
            phscat = true                              # any phantom interaction
        end
        if !is_comp                                    # photoelectric/pair: full absorption
            break
        end
        dir = ndir; E = nE
        if E < egamma_cut                              # below cut: remaining E deposited here
            INSCAN ? (esum += E) : (phscat = true)
            break
        end
    end
    (escaped, pos, dir, E, reached, hx, hy, hz, hiz, hiphi, esum, nblk, pbz, pbphi, phscat)
end

"""
    navigate_single_photons(geom, E0, pos0, dir0, rng; egamma_cut=0.010) -> NamedTuple

The production singles path: transport one photon across the whole geometry exactly as
`navigate_photon` does — same physics, same RNG draw order — but **fold** each interaction
straight into stack-local accumulators instead of building a `Vector{NavStep}`. Returns the
singles summary (`reached, x, y, z, iz, iphi, e, nblocks, phscat`; see `_reduce_steps`) with
**zero heap allocation per photon**, so it scales to the 10⁸-decay production runs where the
vector path's GC would serialise the threads.

It reuses the allocation-free kernels (`sigma_macro`, `sample_distance`, `distance_to_exit`,
`sample_interaction_t`) and `_leaf_reduce` for the per-volume step loop; only the bookkeeping
shell is forked from `propagate_photon` / `navigate_photon`. The tests assert it returns
exactly `_reduce_steps(navigate_photon(...))` on the same seed, guarding the shells against
drift, and that it allocates nothing per photon.
"""
function navigate_single_photons(geom::Geometry, E0::Real, pos0, dir0, rng::AbstractRNG;
                                 egamma_cut::Float64=0.010)
    pos = (Float64(pos0[1]), Float64(pos0[2]), Float64(pos0[3]))
    dir = (Float64(dir0[1]), Float64(dir0[2]), Float64(dir0[3]))
    n   = sqrt(dir[1]^2 + dir[2]^2 + dir[3]^2)
    dir = (dir[1]/n, dir[2]/n, dir[3]/n)
    E   = Float64(E0)

    reached = false
    hx = 0.0; hy = 0.0; hz = 0.0; hiz = -1; hiphi = -1
    esum = 0.0; nblk = 0; pbz = -2; pbphi = -2
    phscat = false

    nseg = 0
    while true
        probe = (pos[1] + NAV_EPS*dir[1], pos[2] + NAV_EPS*dir[2], pos[3] + NAV_EPS*dir[3])
        tag = locate_tag(geom, probe)
        tag === :none && break                         # left the world

        if tag === :world                              # ── air: straight skip, no interaction
            d = next_boundary(geom, geom.world, pos, dir)
            isfinite(d) || break
            pos = (pos[1] + d*dir[1], pos[2] + d*dir[2], pos[3] + d*dir[3])
        elseif tag === :phantom                        # ── phantom leaf
            escaped, pos, dir, E, reached, hx, hy, hz, hiz, hiphi, esum, nblk, pbz, pbphi, phscat =
                _leaf_reduce(geom.phantom, Val(false), nothing, pos, dir, E, egamma_cut,
                             reached, hx, hy, hz, hiz, hiphi, esum, nblk, pbz, pbphi, phscat, rng)
            escaped || break                           # absorbed / below cut → history ends
        else                                           # ── scanner (ring) leaf
            sc = geom.scanner::Scanner
            escaped, pos, dir, E, reached, hx, hy, hz, hiz, hiphi, esum, nblk, pbz, pbphi, phscat =
                _leaf_reduce(sc.volume, Val(true), sc, pos, dir, E, egamma_cut,
                             reached, hx, hy, hz, hiz, hiphi, esum, nblk, pbz, pbphi, phscat, rng)
            escaped || break
        end

        nseg += 1
        nseg > MAX_SEGMENTS && break
    end

    (reached=reached, x=hx, y=hy, z=hz, iz=hiz, iphi=hiphi, e=esum, nblocks=nblk, phscat=phscat)
end
