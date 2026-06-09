#!/usr/bin/env julia
# Step 1 driver: shoot a pencil of 511 keV photons into the phantom and write the
# photon stack (the interactions of each primary) to a CSV.
#
# Run from the repo root:
#   julia --project=. scripts/propagate_gammas_in_phantom.jl --nevents 10000

using PTCryspMC
using ArgParse
using Random

function parse_cli()
    s = ArgParseSettings(description="Propagate pencil 511 keV photons through the phantom; write the photon stack per event.")
    @add_arg_table! s begin
        "--geometry"; help = "geometry JSON"; default = joinpath(@__DIR__, "..", "geometry", "geometry.json")
        "--data";    help = "data dir";      default = joinpath(@__DIR__, "..", "data")
        "--out";     help = "output CSV";    default = joinpath(@__DIR__, "..", "output", "phantom_stack.csv")
        "--nevents"; help = "n photons";     arg_type = Int;     default = 10000
        "--energy";  help = "energy [keV]";  arg_type = Float64; default = 511.0
        "--seed";    help = "RNG seed";      arg_type = Int;     default = 1234
    end
    parse_args(s)
end

function main()
    a = parse_cli()
    mats = load_materials(a["data"])
    geom = load_geometry(a["geometry"], mats)
    pv = geom.phantom
    cyl = solid(pv)

    E0 = a["energy"] / 1000.0           # MeV
    rng = MersenneTwister(a["seed"])

    # Pencil: enter at the centre of the phantom's -z face, travel along +z.
    start = (pv.position[1], pv.position[2], pv.position[3] - cyl.half_length_cm)
    dir   = (0.0, 0.0, 1.0)

    println("phantom '$(name(pv))': $(material(pv).name), cylinder R=$(cyl.radius_cm) cm, " *
            "half-length=$(cyl.half_length_cm) cm; mfp@$(a["energy"]) keV = " *
            "$(round(mfp(material(pv), E0), digits=3)) cm")

    mkpath(dirname(a["out"]))
    nrows = 0
    open(a["out"], "w") do io
        println(io, "event_number,step,x_mm,y_mm,z_mm,e_in_keV,e_dep_keV,process")
        for ev in 1:a["nevents"]
            recs = propagate_photon(E0, start, dir, pv, rng).recs
            for (k, r) in enumerate(recs)
                println(io, join((ev, k,
                    round(r.x * 10, digits=4), round(r.y * 10, digits=4), round(r.z * 10, digits=4),
                    round(r.e_in * 1000, digits=4), round(r.e_dep * 1000, digits=4), r.process), ","))
                nrows += 1
            end
        end
    end
    println("wrote $(a["nevents"]) events ($nrows interaction rows) -> $(a["out"])")
end

main()
