#!/usr/bin/env julia
# QA for a CLINICAL (activity-driven) run: verify the annihilation points in the singles populate
# each source region in proportion to its activity. Rebuilds the ClinicSource from the config, then
# classifies each single's emission point (x0,y0,z0) into the SMALLEST region that contains it
# (an insert wins over the background it sits inside) using the geometry's 3-D `is_inside`, and
# reports, per region, the measured point density against the activity the config specifies. For a
# k:1 hot:background phantom (background filling the inserts too), the density ratio should be ≈ k.
#
#   julia --project=. scripts/tests/check_clinic_regions.jl --config runs/nema_iq_f18_csi.toml

using PTCryspMC
using ArgParse
using Printf

function parse_cli()
    s = ArgParseSettings(description="Check a clinical run's per-region draw against its activity.")
    @add_arg_table! s begin
        "--config";  help = "run config TOML (must be [source].mode=clinic)"; required = true
        "--singles"; help = "override the singles path (default prod/<tag>/singles.h5)"; default = ""
    end
    parse_args(s)
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    String(cfg_get(cfg, "source", "mode", "")) == "clinic" ||
        error("[source].mode must be 'clinic' for this check")
    tag    = run_tag(cfg, a["config"])
    outdir = joinpath(rp(prod_base(cfg)), tag)
    singles = isempty(a["singles"]) ? joinpath(outdir, "singles.h5") : rp(a["singles"])
    isfile(singles) || error("singles file '$singles' not found")

    mats = load_materials(rp(cfg_get(cfg, "paths", "data", "data")))
    geom = load_geometry(rp(cfg_get(cfg, "geometry", "file", "geometry/geometry.json")), mats)
    cs   = load_clinic_source(cfg, geom)

    n      = length(cs.regions)
    vols   = [volume(pv) for pv in cs.regions]                # cm³ (= mL)
    act    = [cs.conc[k] * vols[k] for k in 1:n]              # Bq per region
    order  = sortperm(vols)                                   # smallest first → inserts before background
    counts = zeros(Int, n)
    tot    = Ref(0)

    foreach_singles_hdf5(singles) do b
        for i in 1:length(b)
            p = (decode_xyz(b.x0[i]) / 10, decode_xyz(b.y0[i]) / 10, decode_xyz(b.z0[i]) / 10)  # mm → cm
            tot[] += 1
            for k in order
                if is_inside(cs.regions[k], p)
                    counts[k] += 1
                    break
                end
            end
        end
    end

    big   = order[end]                                        # the largest region = background
    dens  = [counts[k] / vols[k] for k in 1:n]                # measured point density [1/mL]
    densb = dens[big]
    @printf("run '%s': %d singles over %d region(s); total activity %.1f kBq\n",
            tag, tot[], n, cs.A0 / 1e3)
    println("  region        conc[kBq/mL]   vol[mL]    count    density/bg")
    for k in order
        @printf("  %-12s  %8.2f     %8.1f  %9d   %7.2f\n",
                name(cs.regions[k]), cs.conc[k] / 1e3, vols[k], counts[k],
                densb > 0 ? dens[k] / densb : NaN)
    end
    println("  (density/bg ≈ effective hot:background ratio; background row = 1.00 by construction)")
end

main()
