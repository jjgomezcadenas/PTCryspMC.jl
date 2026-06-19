#!/usr/bin/env julia
# Why does the water phantom drop the clean-coincidence fraction so much? This script
# isolates the phantom's effect by running the SAME back-to-back source two ways and
# comparing, per crystal material:
#
#   air-only  — each photon flies straight through air to the ring (the unit-test path,
#               shoot_into_ring.jl): no phantom, full 511 keV arrives.
#   phantom   — each photon is NAVIGATED from the phantom centre through water → air →
#               ring (navigate_photon): ~8 cm of water in the way.
#
# The metric is the BOTH-PHOTON FULL-ENERGY fraction (truth): both photons of an
# annihilation deposit ≥505 keV in a single crystal each. This is a TRUTH-level cut (no
# detector resolution): at truth a fully-absorbed photon deposits exactly 511 keV, so it
# isolates the full-energy peak — but a ≥505 keV cut also rejects every phantom-scattered
# photon, so this fraction is the UNSCATTERED, fully-contained subset. It is NOT the
# coincidence efficiency: the real selection (Step 4) smears by FWHM(E) and applies a
# ±2·FWHM energy window that also keeps scattered photons as *scatter* coincidences.
#
# Reported under two conditionings so the air↔phantom comparison is apples-to-apples:
#   • over ALL events           — the both-photon full-energy efficiency;
#   • among BOTH-REACHED events — removes the geometric acceptance, isolating the
#                                 energy/containment loss.
#
# The expected decomposition (and what the numbers confirm):
#   a 511 keV photon survives 8 cm of water unscattered with probability
#   exp(−8/mfp) ≈ exp(−8/10.44) ≈ 0.46; a photon that Compton-scatters in the phantom
#   loses energy and falls below the 505 keV cut. Hence
#       fullE_phantom/γ ≈ P(unscattered) × fullE_air/γ,
#   and since back-to-back photons share |cosθ| (acceptance is 100% correlated) but
#   scatter independently,
#       fullE_phantom(both-reached) ≈ P(unscattered)² × fullE_air(both-reached).
#
# Run from the repo root:
#   julia --project=. scripts/phantom_effect_on_coincidences.jl --nevents 30000
#   julia --project=. scripts/phantom_effect_on_coincidences.jl --materials CsI,BGO

using PTCryspMC
using ArgParse
using Random
using Printf

function parse_cli()
    s = ArgParseSettings(description="Compare clean-coincidence fraction air-only vs through the water phantom.")
    @add_arg_table! s begin
        "--data";      help = "data dir";        default = joinpath(@__DIR__, "..", "..", "data")
        "--geometry";  help = "geometry JSON";   default = joinpath(@__DIR__, "..", "..", "geometry", "geometry.json")
        "--materials"; help = "comma-separated crystal materials"; default = "CsI,BGO"
        "--nevents";   help = "n annihilations"; arg_type = Int;     default = 30000
        "--cutoff";    help = "low-energy cutoff [keV]"; arg_type = Float64; default = 10.0
        "--seed";      help = "RNG seed";        arg_type = Int;     default = 1234
    end
    parse_args(s)
end

"A random unit vector, isotropic on the sphere."
@inline function rand_unit(rng)
    c = 2.0 * rand(rng) - 1.0
    ϕ = 2π * rand(rng)
    s = sqrt(max(0.0, 1.0 - c^2))
    (s * cos(ϕ), s * sin(ϕ), c)
end

"Containment of one photon's scanner deposits: full 511 keV in a single crystal."
function scanner_outcome(Escan::Float64, nblk::Int)
    reached   = nblk >= 1
    contained = Escan >= 0.505 && nblk == 1
    (reached, contained)
end

# One photon, air-only: straight through air to the ring, then transport in the crystal.
function gamma_air(sc::Scanner, dir, E0, cut, rng)
    de = distance_to_entry((0.0, 0.0, 0.0), dir, sc.volume)
    isfinite(de) || return (false, false, false)        # missed the ring
    entry = (de*dir[1], de*dir[2], de*dir[3])
    recs = propagate_photon(E0, entry, dir, sc.volume, rng; egamma_cut=cut).recs
    Escan = 0.0; blocks = Set{Tuple{Int,Int}}()
    for r in recs
        r.process == :escape && continue
        Escan += r.e_dep
        push!(blocks, block_index(sc, (r.x, r.y, r.z)))
    end
    (scanner_outcome(Escan, length(blocks))..., false)   # no phantom: never scattered there
end

# One photon, navigated from the phantom centre through water → air → ring.
function gamma_phantom(geom::Geometry, dir, E0, cut, rng)
    steps = navigate_photon(geom, E0, (0.0, 0.0, 0.0), dir, rng; egamma_cut=cut)
    Escan = 0.0; blocks = Set{Tuple{Int,Int}}(); phscat = false
    for st in steps
        st.hit.process == :escape && continue
        if st.volume == :phantom
            phscat = true
        elseif st.volume == :scanner
            Escan += st.hit.e_dep
            push!(blocks, (st.iz, st.iphi))
        end
    end
    (scanner_outcome(Escan, length(blocks))..., phscat)
end

# Run n annihilations one way; `emit` directions come from `rng_emit` (reset identically
# per mode, so both modes see the SAME annihilations), transport from `rng_trk`.
function tally(gamma_fn, n, E0, cut, rng_emit, rng_trk)
    nreach = 0; nfe = 0; nunscat = 0
    nboth_reach = 0; both_fe_all = 0; both_fe_reach = 0
    for _ in 1:n
        d = rand_unit(rng_emit)
        r = (false, false); fe = (false, false)
        for (i, dir) in enumerate((d, (-d[1], -d[2], -d[3])))
            reached, full_E, phscat = gamma_fn(dir, E0, cut, rng_trk)
            nreach += reached; nfe += full_E; nunscat += !phscat
            r  = i == 1 ? (reached, r[2]) : (r[1], reached)
            fe = i == 1 ? (full_E, fe[2]) : (fe[1], full_E)
        end
        both_fe_all += (fe[1] && fe[2])
        if r[1] && r[2]
            nboth_reach += 1
            both_fe_reach += (fe[1] && fe[2])
        end
    end
    ng = 2n
    (reach=nreach/ng, fe_all=nfe/ng, fe_reach=nfe/max(nreach, 1), unscat=nunscat/ng,
     bothfe_all=both_fe_all/n, bothfe_reach=both_fe_reach/max(nboth_reach, 1))
end

function main()
    a = parse_cli()
    mats = load_materials(a["data"])
    g0   = load_geometry(a["geometry"], mats)
    g0.scanner === nothing && error("geometry $(a["geometry"]) has no scanner section")
    E0  = 0.511
    cut = a["cutoff"] / 1000.0
    n   = a["nevents"]

    pmat = material(g0.phantom)
    mfp_w = mfp(pmat, E0)
    Rp = solid(g0.phantom).radius_cm
    psurv = exp(-Rp / mfp_w)
    @printf("phantom: %s, R=%.1f cm, mfp@511keV=%.2f cm → P(survive %.0f cm unscattered)=%.1f%%\n",
            name(g0.phantom), Rp, mfp_w, Rp, 100psurv)
    println("back-to-back 511 keV from the phantom centre; $n events per material\n")

    println("metric = both-photon FULL-ENERGY (truth, ≥505 keV in 1 crystal each) = the")
    println("UNSCATTERED, fully-contained subset — NOT the coincidence efficiency.\n")
    hdr = ("material", "mode", "reach/γ", "fullE/γ", "fullE/γ|R", "unscat/γ", "both-fE(all)", "both-fE(R)")
    @printf("%-8s %-8s %8s %8s %9s %9s %13s %11s\n", hdr...)
    println("-"^82)

    for m in split(a["materials"], ',')
        m = strip(m); haskey(mats, m) || error("material '$m' not found")
        sc0 = g0.scanner
        sc  = Scanner(PhysicalVolume(LogicalVolume(name(sc0.volume), solid(sc0.volume),
                                                   mats[m]), sc0.volume.position),
                      sc0.n_phi, sc0.n_z)
        geom = Geometry(g0.world, g0.phantom, sc)

        air = tally((dir, E, c, rng) -> gamma_air(sc, dir, E, c, rng), n, E0, cut,
                    MersenneTwister(a["seed"]), MersenneTwister(a["seed"] + 1))
        phn = tally((dir, E, c, rng) -> gamma_phantom(geom, dir, E, c, rng), n, E0, cut,
                    MersenneTwister(a["seed"]), MersenneTwister(a["seed"] + 1))

        for (lbl, r) in (("air-only", air), ("phantom", phn))
            @printf("%-8s %-8s %7.1f%% %7.1f%% %8.1f%% %8.1f%% %12.1f%% %10.1f%%\n",
                    m, lbl, 100r.reach, 100r.fe_all, 100r.fe_reach, 100r.unscat,
                    100r.bothfe_all, 100r.bothfe_reach)
        end
        # Check the decomposition: both-fE_phantom(R) ≈ unscat² × both-fE_air(R).
        pred = phn.unscat^2 * air.bothfe_reach
        @printf("  → expect both-fE(R)_phantom ≈ unscat²·both-fE(R)_air = %.2f²·%.1f%% = %.1f%%  (measured %.1f%%)\n\n",
                phn.unscat, 100air.bothfe_reach, 100pred, 100phn.bothfe_reach)
    end
end

main()
