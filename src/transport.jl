# Photon-only transport through a single cylinder, recording the interaction stack.
# We follow only the photon: at a Compton interaction the recoil energy is deposited
# locally and the photon continues; photoelectric deposits the full energy and ends
# the history. (See docs/PTCryspMC_phys.tex.)

"One recorded interaction of a photon history. Positions [cm], energies [MeV]."
struct Interaction
    x::Float64
    y::Float64
    z::Float64
    e_in::Float64      # photon energy entering this interaction
    e_dep::Float64     # energy deposited here
    process::Symbol    # :compton, :photoelectric, :pair, :below_cut, :escape
end

"""
    Transported

The result of transporting a photon through one volume: the interaction stack `recs`
plus the photon's exit state — `pos`, `dir`, surviving energy `E`, and `escaped`
(`true` if it left the volume alive, `false` if it was absorbed or fell below cut).
The exit state is what lets `navigate_photon` carry the photon into the next volume.
"""
struct Transported
    recs::Vector{Interaction}
    pos::NTuple{3,Float64}
    dir::NTuple{3,Float64}
    E::Float64
    escaped::Bool
end

"""
    propagate_photon(E0_MeV, pos0, dir0, pv, rng; egamma_cut=0.010) -> Transported

Transport one photon of energy `E0_MeV` from `pos0` along `dir0` through the single
physical volume `pv` (its solid filled with its material). Returns a `Transported`:
the stack of interactions — Compton scatters (recoil deposited, photon continues), a
terminating photoelectric absorption, an `:escape` record at the exit point when the
photon leaves the volume, or a `:below_cut` record when a scattered photon falls below
`egamma_cut` [MeV] — together with the exit state. The interaction physics is shared
with `navigate_photon` through `sample_interaction`. Callers wanting only the stack
read `.recs`.
"""
function propagate_photon(E0_MeV::Real, pos0, dir0, pv::PhysicalVolume,
                          rng::AbstractRNG; egamma_cut::Float64=0.010)::Transported
    E = Float64(E0_MeV)
    pos = collect(Float64, pos0)
    dir = collect(Float64, dir0)
    dir ./= sqrt(sum(abs2, dir))
    mat = material(pv)

    recs = Interaction[]
    while true
        ΣC, ΣPh, ΣP = sigma_macro(mat, E)
        Σ = ΣC + ΣPh + ΣP
        s = Σ > 0.0 ? sample_distance(Σ, rng) : Inf
        d_exit = distance_to_exit(pos, dir, pv)

        if s >= d_exit
            pos = pos .+ d_exit .* dir
            push!(recs, Interaction(pos[1], pos[2], pos[3], E, 0.0, :escape))
            return Transported(recs, (pos[1], pos[2], pos[3]), (dir[1], dir[2], dir[3]), E, true)
        end

        pos = pos .+ s .* dir
        proc, e_dep, ndir, nE = sample_interaction(E, dir, ΣC, ΣPh, ΣP, rng)
        push!(recs, Interaction(pos[1], pos[2], pos[3], E, e_dep, proc))

        if proc !== :compton
            # photoelectric (or pair, negligible at 511 keV): full absorption, history ends
            return Transported(recs, (pos[1], pos[2], pos[3]), (dir[1], dir[2], dir[3]), 0.0, false)
        end

        dir = ndir
        E = nE
        if E < egamma_cut
            push!(recs, Interaction(pos[1], pos[2], pos[3], E, E, :below_cut))
            return Transported(recs, (pos[1], pos[2], pos[3]), (dir[1], dir[2], dir[3]), 0.0, false)
        end
    end
end
