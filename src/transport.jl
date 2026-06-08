# Photon-only transport through a single cylinder, recording the interaction stack.
# We follow only the photon: at a Compton interaction the recoil energy is deposited
# locally and the photon continues; photoelectric deposits the full energy and ends
# the history. (See docs/pet_simulation.tex.)

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
    propagate_photon(E0_MeV, pos0, dir0, pv, rng; egamma_cut=0.010)
        -> Vector{Interaction}

Transport one photon of energy `E0_MeV` from `pos0` along `dir0` through the
physical volume `pv` (its solid filled with its material). Returns the stack of
interactions: Compton scatters (recoil deposited, photon continues), a
terminating photoelectric absorption, an `:escape` record at the exit point when
the photon leaves the volume, or a `:below_cut` record when a scattered photon
falls below `egamma_cut` [MeV].
"""
function propagate_photon(E0_MeV::Real, pos0, dir0, pv::PhysicalVolume,
                          rng::AbstractRNG; egamma_cut::Float64=0.010)::Vector{Interaction}
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
            ep = pos .+ d_exit .* dir
            push!(recs, Interaction(ep[1], ep[2], ep[3], E, 0.0, :escape))
            break
        end

        pos = pos .+ s .* dir
        proc = sample_process(ΣC / Σ, ΣPh / Σ, ΣP / Σ, rng)

        if proc === :compton
            Eprime, cosθ = sample_compton(E, rng)
            push!(recs, Interaction(pos[1], pos[2], pos[3], E, E - Eprime, :compton))
            ϕ = 2π * rand(rng)
            sinθ = sqrt(max(0.0, 1.0 - cosθ^2))
            dir = rotate_to_global(Float64[sinθ*cos(ϕ), sinθ*sin(ϕ), cosθ], dir)
            E = Eprime
            if E < egamma_cut
                push!(recs, Interaction(pos[1], pos[2], pos[3], E, E, :below_cut))
                break
            end
        else
            # photoelectric (or pair, negligible at 511 keV): full absorption, history ends
            push!(recs, Interaction(pos[1], pos[2], pos[3], E, E, proc))
            break
        end
    end
    recs
end
