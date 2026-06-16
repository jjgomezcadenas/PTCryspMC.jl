#!/usr/bin/env julia
# Per-event timing of the back-to-back transport (the work in
# shoot_into_ring.jl, minus the CSV writing): emit an isotropic
# 511 keV pair from the origin, navigate each photon straight to the ring, transport
# it through the crystal, and tag every interaction with its block. Reports the
# compute cost per annihilation, after a warm-up, for each material.
#
# Run from the repo root:
#   julia --project=. scripts/bench_back_to_back.jl --nevents 500000 --materials CsI,BGO

using PTCryspMC
using ArgParse
using Random

function parse_cli()
    s = ArgParseSettings(description="Time the back-to-back 511 keV transport per annihilation.")
    @add_arg_table! s begin
        "--data";      help = "data dir";        default = joinpath(@__DIR__, "..", "data")
        "--geometry";  help = "geometry JSON";   default = joinpath(@__DIR__, "..", "geometry", "geometry.json")
        "--materials"; help = "comma-separated crystal materials"; default = "CsI,BGO"
        "--nevents";   help = "n annihilations"; arg_type = Int;     default = 500000
        "--cutoff";    help = "low-energy cutoff [keV]"; arg_type = Float64; default = 10.0
        "--seed";      help = "RNG seed";        arg_type = Int;     default = 1234
    end
    parse_args(s)
end

@inline function rand_unit(rng)
    c = 2.0 * rand(rng) - 1.0
    ϕ = 2π * rand(rng)
    s = sqrt(max(0.0, 1.0 - c^2))
    (s * cos(ϕ), s * sin(ϕ), c)
end

# The per-event core, no IO. Returns the photon hit count so the work isn't elided.
function run!(sc::Scanner, n::Int, cut_MeV::Float64, rng)
    E0 = 0.511
    acc = 0
    @inbounds for _ in 1:n
        d = rand_unit(rng)
        for dir in (d, (-d[1], -d[2], -d[3]))
            de = distance_to_entry((0.0, 0.0, 0.0), dir, sc.volume)
            isfinite(de) || continue
            acc += 1
            entry = (de * dir[1], de * dir[2], de * dir[3])
            recs = propagate_photon(E0, entry, dir, sc.volume, rng; egamma_cut=cut_MeV).recs
            for r in recs
                iϕ, iz = block_index(sc, (r.x, r.y, r.z))
                acc += (iϕ + iz) & 0   # touch the indices without changing the count
            end
        end
    end
    acc
end

function main()
    a = parse_cli()
    mats = load_materials(a["data"])
    geom = load_geometry(a["geometry"], mats)
    geom.scanner === nothing && error("geometry $(a["geometry"]) has no scanner section")
    sc0 = geom.scanner
    cut = a["cutoff"] / 1000.0
    n   = a["nevents"]

    println("back-to-back 511 keV, cutoff $(a["cutoff"]) keV, ring $(name(sc0.volume)) " *
            "($(nblocks(sc0)) blocks); $n events per material")
    for m in split(a["materials"], ',')
        m = strip(m)
        haskey(mats, m) || error("material '$m' not found")
        sc  = Scanner(PhysicalVolume(LogicalVolume(name(sc0.volume), solid(sc0.volume), mats[m]),
                                     sc0.volume.position), sc0.n_phi, sc0.n_z)
        run!(sc, min(n, 1000), cut, MersenneTwister(a["seed"]))   # warm up / compile
        t = @elapsed run!(sc, n, cut, MersenneTwister(a["seed"]))
        println("  $(rpad(m, 5)): $(round(t, digits=3)) s  =  " *
                "$(round(t/n*1e6, digits=3)) µs/event  " *
                "($(round(n/t/1e6, digits=2)) M events/s)")
    end
end

main()
