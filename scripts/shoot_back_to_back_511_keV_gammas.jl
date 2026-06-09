#!/usr/bin/env julia
# Unit test: a point source at the origin emits back-to-back 511 keV photon pairs
# isotropically into air; each photon flies straight to the detector ring (or
# misses it) and is transported through the crystal. We record the full interaction
# stack of both photons, tagged with the crystal (iz, iphi) of each interaction.
#
# Emission is script-local (no phantom, no scenario): a 1 mm "test sphere" of air is
# geometrically inert, so the source is just a point (optional --source-radius). The
# water phantom in geometry.json is not used — the interior is treated as air, so
# every photon reaches the ring or escapes with no scattering on the way.
#
# Run from the repo root:
#   julia --project=. scripts/shoot_back_to_back_511_keV_gammas.jl --nevents 100000

using PTCryspMC
using ArgParse
using Random

function parse_cli()
    s = ArgParseSettings(description="Back-to-back 511 keV pairs from a point source into the detector ring; write both photon stacks.")
    @add_arg_table! s begin
        "--data";          help = "data dir";        default = joinpath(@__DIR__, "..", "data")
        "--geometry";      help = "geometry JSON";   default = joinpath(@__DIR__, "..", "geometry", "geometry.json")
        "--out";           help = "output CSV (default output/b2b_<material>_stack.csv)"; default = ""
        "--material";      help = "crystal material"; default = "CsI"
        "--nevents";       help = "n annihilations"; arg_type = Int;     default = 100000
        "--energy";        help = "energy [keV]";    arg_type = Float64; default = 511.0
        "--cutoff";        help = "low-energy cutoff [keV] (default 10 = XCOM min)"; arg_type = Float64; default = 10.0
        "--source-radius"; help = "source sphere radius [mm] (0 = point)"; arg_type = Float64; default = 0.0
        "--seed";          help = "RNG seed";        arg_type = Int;     default = 1234
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

"""
    emit_pair(rng, radius_cm) -> (pos, dir1, dir2)

One annihilation: a position (the origin, or uniform within a sphere of `radius_cm`)
and two exactly back-to-back isotropic directions. The seed of a future `Source`.
"""
function emit_pair(rng, radius_cm::Float64)
    d = rand_unit(rng)
    pos = if radius_cm > 0.0
        r = radius_cm * cbrt(rand(rng))            # uniform in the sphere
        u = rand_unit(rng)
        (r * u[1], r * u[2], r * u[3])
    else
        (0.0, 0.0, 0.0)
    end
    (pos, d, (-d[1], -d[2], -d[3]))
end

function main()
    a = parse_cli()
    mats = load_materials(a["data"])
    geom = load_geometry(a["geometry"], mats)
    geom.scanner === nothing && error("geometry $(a["geometry"]) has no scanner section")

    # Take the ring geometry/segmentation from the scanner, but use the chosen material.
    haskey(mats, a["material"]) || error("material '$(a["material"])' not found")
    sc0 = geom.scanner
    cs  = solid(sc0.volume)
    sc  = Scanner(PhysicalVolume(LogicalVolume(name(sc0.volume), cs, mats[a["material"]]),
                                 sc0.volume.position), sc0.n_phi, sc0.n_z)

    out = isempty(a["out"]) ?
        joinpath(@__DIR__, "..", "output", "b2b_$(lowercase(a["material"]))_stack.csv") : a["out"]
    E0       = a["energy"] / 1000.0          # MeV
    cut_MeV  = a["cutoff"] / 1000.0
    radius   = a["source-radius"] / 10.0     # mm -> cm
    rng      = MersenneTwister(a["seed"])

    println("source: back-to-back $(a["energy"]) keV, isotropic, radius $(a["source-radius"]) mm; " *
            "ring $(name(sc.volume)) ($(a["material"])), nblocks=$(nblocks(sc)), cutoff $(a["cutoff"]) keV")

    n_hit = 0; n_miss = 0
    mkpath(dirname(out))
    nrows = 0
    open(out, "w") do io
        println(io, "event_number,gamma,step,x_mm,y_mm,z_mm,e_in_keV,e_dep_keV,process,iz,iphi")
        for ev in 1:a["nevents"]
            pos, d1, d2 = emit_pair(rng, radius)
            for (g, dir) in ((1, d1), (2, d2))
                de = distance_to_entry(pos, dir, sc.volume)   # straight through air to the ring
                if !isfinite(de)
                    n_miss += 1
                    continue
                end
                n_hit += 1
                entry = (pos[1] + de*dir[1], pos[2] + de*dir[2], pos[3] + de*dir[3])
                recs = propagate_photon(E0, entry, dir, sc.volume, rng; egamma_cut=cut_MeV).recs
                for (k, r) in enumerate(recs)
                    iϕ, iz = block_index(sc, (r.x, r.y, r.z))
                    println(io, join((ev, g, k,
                        round(r.x*10, digits=4), round(r.y*10, digits=4), round(r.z*10, digits=4),
                        round(r.e_in*1000, digits=4), round(r.e_dep*1000, digits=4), r.process, iz, iϕ), ","))
                    nrows += 1
                end
            end
        end
    end

    ngamma = 2 * a["nevents"]
    println("wrote $(a["nevents"]) events ($nrows interaction rows) -> $out")
    println("ring acceptance: $(round(100*n_hit/ngamma, digits=1))% of photons hit " *
            "($n_hit hits, $n_miss misses of $ngamma)")
end

main()
