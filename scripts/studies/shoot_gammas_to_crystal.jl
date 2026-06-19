#!/usr/bin/env julia
# Shoot 511 keV photons into a single scintillator crystal (a box) and write the
# interaction stack per event, for a containment / topology study. The crystal
# material is set by --material (e.g. CsI, BGO).
#
# Geometry: a crystal box of width x width x depth mm, placed so its entry face is
# at world z = 0 and the crystal spans z in [0, depth]. So an interaction's z is
# its depth from the entry face.
#
# Source (fully parametric):
#   --beam-xy <mm>      half-width of the uniform start square on the entry face;
#                       0 = central pencil, width/2 = the whole face.
#   --beam-opening <deg> TOTAL angular opening about +z; the photon's polar angle
#                       from +z is uniform in [0, opening/2] with random azimuth
#                       (a 3-D cone of half-angle opening/2). 0 = parallel to z.
#   --tag <name>        output suffix, so runs never clobber:
#                       output/<material>_crystal_<tag>_stack.csv
#
# Run from the repo root:
#   julia --project=. scripts/shoot_gammas_to_crystal.jl --material CsI --tag cone   --beam-xy 24 --beam-opening 45
#   julia --project=. scripts/shoot_gammas_to_crystal.jl --material BGO --tag pencil --beam-xy 0  --beam-opening 0

using PTCryspMC
using ArgParse
using Random

function parse_cli()
    s = ArgParseSettings(description="Shoot 511 keV photons into a CsI crystal box; write the interaction stack.")
    @add_arg_table! s begin
        "--data";         help = "data dir";          default = joinpath(@__DIR__, "..", "..", "data")
        "--tag";          help = "output suffix (csi_crystal_<tag>_stack.csv)"; default = "cone"
        "--out";          help = "output CSV (overrides --tag)"; default = ""
        "--material";     help = "crystal material";  default = "CsI"
        "--width";        help = "face width [mm]";    arg_type = Float64; default = 48.0
        "--depth";        help = "crystal depth [mm]"; arg_type = Float64; default = 37.0
        "--beam-xy";      help = "start half-width on the face [mm]; 0 = centre"; arg_type = Float64; default = 24.0
        "--beam-opening"; help = "total angular opening about +z [deg]"; arg_type = Float64; default = 45.0
        "--nevents";      help = "n photons";          arg_type = Int;     default = 100000
        "--energy";       help = "energy [keV]";       arg_type = Float64; default = 511.0
        "--seed";         help = "RNG seed";           arg_type = Int;     default = 1234
    end
    parse_args(s)
end

function main()
    a = parse_cli()
    mat = load_material(a["data"], a["material"])
    out = isempty(a["out"]) ?
        joinpath(@__DIR__, "..", "..", "output",
                 "$(lowercase(a["material"]))_crystal_$(a["tag"])_stack.csv") : a["out"]

    # Box centred at the origin, placed so the entry (-z) face sits at world z = 0.
    hx = a["width"] / 2 / 10.0     # mm -> cm half-widths
    hy = hx
    hz = a["depth"] / 2 / 10.0
    box = Box(hx, hy, hz)
    pv = PhysicalVolume(LogicalVolume(a["material"], box, mat), (0.0, 0.0, hz))

    bxy = a["beam-xy"] / 10.0                  # mm -> cm start half-width
    bxy <= hx || error("--beam-xy ($(a["beam-xy"]) mm) exceeds the face half-width ($(a["width"]/2) mm)")
    E0 = a["energy"] / 1000.0                  # MeV
    θhalf = deg2rad(a["beam-opening"] / 2)     # cone half-angle = opening / 2
    rng = MersenneTwister(a["seed"])

    println("crystal '$(a["material"])': $(a["width"])x$(a["width"])x$(a["depth"]) mm; " *
            "$(a["energy"]) keV; beam-xy ±$(a["beam-xy"]/2) mm, opening $(a["beam-opening"])° " *
            "(half-angle $(round(rad2deg(θhalf), digits=1))°); att.length = $(round(mfp(mat, E0), digits=3)) cm")

    n_contained = 0; n_photo = 0; n_compton = 0
    mkpath(dirname(out))
    nrows = 0
    open(out, "w") do io
        println(io, "event_number,step,x_mm,y_mm,z_mm,e_in_keV,e_dep_keV,process")
        for ev in 1:a["nevents"]
            # entry point on the z = 0 face, uniform over the ±bxy square (bxy=0 → centre)
            x = (2.0 * rand(rng) - 1.0) * bxy
            y = (2.0 * rand(rng) - 1.0) * bxy
            start = (x, y, 0.0)
            # direction in a 3-D cone of half-angle θhalf about +z (always into the crystal)
            θ = θhalf * rand(rng)
            ϕ = 2π * rand(rng)
            dir = (sin(θ) * cos(ϕ), sin(θ) * sin(ϕ), cos(θ))

            recs = propagate_photon(E0, start, dir, pv, rng).recs
            for (k, r) in enumerate(recs)
                println(io, join((ev, k,
                    round(r.x * 10, digits=4), round(r.y * 10, digits=4), round(r.z * 10, digits=4),
                    round(r.e_in * 1000, digits=4), round(r.e_dep * 1000, digits=4), r.process), ","))
                nrows += 1
            end

            # running summary: contained = photon fully stopped (did not escape)
            if recs[end].process !== :escape
                n_contained += 1
                nsc = count(r -> r.process === :compton, recs)
                nsc == 0 ? (n_photo += 1) : (n_compton += 1)
            end
        end
    end

    N = a["nevents"]
    println("wrote $N events ($nrows interaction rows) -> $out")
    println("contained = $(round(100*n_contained/N, digits=1))%  " *
            "(of contained: photo $(round(100*n_photo/max(n_contained,1), digits=1))%, " *
            "compton $(round(100*n_compton/max(n_contained,1), digits=1))%)")
end

main()
