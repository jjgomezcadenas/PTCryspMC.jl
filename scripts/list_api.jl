#!/usr/bin/env julia
# Inventory of the API (Proton Activity) source: the run configs (runs/*.toml with
# [source].mode="api"), the frozen ptcryspg4 scenarios they reference, the produced prod/<tag>/
# outputs, and the repo machinery that implements it. A "where is every part" manifest — read-only.
#
#   julia --project=. scripts/list_api.jl
#   julia --project=. scripts/list_api.jl --scenarios-root ~/Projects/ptcrysp-scenarios/scenarios

using PTCryspMC
using ArgParse
using Printf

function parse_cli()
    s = ArgParseSettings(description="List the API run configs, scenarios, outputs and components.")
    @add_arg_table! s begin
        "--runs";            help = "directory of run configs"; default = "runs"
        "--scenarios-root";  help = "also list every frozen scenario under this dir"; default = ""
    end
    parse_args(s)
end

human(bytes) = bytes < 1e3 ? @sprintf("%d B", bytes) :
               bytes < 1e6 ? @sprintf("%.1f KB", bytes/1e3) :
               bytes < 1e9 ? @sprintf("%.1f MB", bytes/1e6) : @sprintf("%.2f GB", bytes/1e9)

# One data row of a small scenario CSV → a name=>value map (nothing if the file is absent).
function meta_row(path)
    isfile(path) || return nothing
    lines = filter(!isempty, strip.(readlines(path)))
    length(lines) >= 2 || return nothing
    h = strip.(split(lines[1], ',')); v = strip.(split(lines[2], ','))
    Dict(String(h[i]) => String(v[i]) for i in eachindex(h))
end

function scenario_summary(dir, rp)
    d = rp(dir)
    isdir(d) || return "  scenario: $dir  [NOT FOUND on this machine]"
    rm = meta_row(joinpath(d, "run_meta.csv"))
    isos = String[]
    if isfile(joinpath(d, "isotopes.csv"))
        for l in Iterators.drop(eachline(joinpath(d, "isotopes.csv")), 1)
            isempty(strip(l)) && continue
            push!(isos, String(strip(split(l, ',')[2])))
        end
    end
    budgets = sort([match(r"sampling_budget_(.+)\.csv", basename(f))[1]
                    for f in readdir(d; join=true)
                    if occursin(r"sampling_budget_[^_]+\.csv$", basename(f))])
    nemit = isfile(joinpath(d, "emitters.csv")) ? countlines(joinpath(d, "emitters.csv")) - 1 : -1
    io = IOBuffer()
    println(io, "  scenario: $(basename(rstrip(d,'/')))  [FOUND]")
    rm !== nothing && println(io, "    geometry $(get(rm,"geometry","?")), phantom $(get(rm,"phantom_material","?")), ",
                             "n_protons $(get(rm,"n_protons","?"))")
    println(io, "    budgets: $(join(budgets, " "))   isotopes: $(join(isos, " "))   emitters: $nemit")
    String(take!(io))
end

function prod_status(outdir)
    isdir(outdir) || return "    prod: (not produced)"
    io = IOBuffer(); println(io, "    prod/$(basename(outdir))/:")
    for (f, label) in (("singles.h5","nevents"), ("lors_truth.h5","nrows"),
                       ("randoms.h5","nrows"), ("lors_det.h5","nrows"))
        p = joinpath(outdir, f)
        if isfile(p)
            n = Int(singles_hdf5_attr(p, label, -1))
            tag = f == "lors_det.h5" ? "  <- DELIVERABLE" : ""
            @printf(io, "      %-14s %9s   %s %s%s\n", f, human(filesize(p)),
                    label, n >= 0 ? string(n) : "?", tag)
        end
    end
    String(take!(io))
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfgs = sort(filter(f -> endswith(f, ".toml"), readdir(rp(a["runs"]); join=true)))
    apis = Tuple{String,Any}[]
    for f in cfgs
        cfg = read_config(f)
        String(cfg_get(cfg, "source", "mode", "")) == "api" && push!(apis, (f, cfg))
    end

    println("="^78)
    println("API run configs ($(rp(a["runs"]))/*.toml with [source].mode = \"api\"):  $(length(apis)) found")
    println("="^78)
    scen_dirs = String[]
    for (f, cfg) in apis
        tag  = run_tag(cfg, f)
        scn  = String(cfg_get(cfg, "source", "scenario_dir", ""));  push!(scen_dirs, scn)
        found = isdir(rp(scn)) ? "" : "  [scenario NOT FOUND]"
        println("\n▸ $tag   ($(basename(f)))")
        println("    detector : $(cfg_get(cfg,"transport","crystal_material","CsI"))   " *
                "geometry: $(cfg_get(cfg,"geometry","file","?"))")
        println("    scenario : $(basename(rstrip(scn,'/')))  ($scn)$found")
        @printf("    source   : budget=%s dose=%sGy master_seed=%s realization=%s keep_escaped=%s\n",
                cfg_get(cfg,"source","budget","fast"), cfg_get(cfg,"source","dose_Gy",1.0),
                cfg_get(cfg,"source","master_seed",1), cfg_get(cfg,"source","realization",0),
                cfg_get(cfg,"source","keep_escaped",false))
        @printf("    response : sigma_xyz=%smm eres=%s reco_emin=%skeV tau=%sns\n",
                cfg_get(cfg,"detector","sigma_xyz_mm",0.0), cfg_get(cfg,"detector","eres","(crystal)"),
                cfg_get(cfg,"detector","reco_emin_keV",cfg_get(cfg,"detector","emin_keV",0.0)),
                cfg_get(cfg,"timing","tau_ns",0.0))
        print(prod_status(joinpath(rp(prod_base(cfg)), tag)))
    end

    println("\n" * "="^78); println("Referenced scenarios:"); println("="^78)
    for scn in unique(scen_dirs); print(scenario_summary(scn, rp)); end
    if !isempty(a["scenarios-root"]) && isdir(a["scenarios-root"])
        println("\nAll scenarios under $(a["scenarios-root"]):")
        for d in sort(readdir(a["scenarios-root"]; join=true))
            isdir(d) && isfile(joinpath(d, "emitters.csv")) && print(scenario_summary(d, identity))
        end
    end

    println("\n" * "="^78); println("Components (the API machinery, tracked in the repo):"); println("="^78)
    comps = [
        ("engine",   ["src/scenario.jl", "src/geometry.jl", "src/source.jl", "src/singles.jl"]),
        ("data",     ["data/materials.json", "data/xcom_brain.csv"]),
        ("geometry", ["geometry/geometry_head.json"]),
        ("chain",    ["scripts/simulate_source_mt.jl", "scripts/build_true_coincidences_from_singles.jl",
                      "scripts/build_randoms_from_singles.jl", "scripts/reco_lors.jl", "scripts/run/run_prod.sh"]),
        ("QA",       ["scripts/tests/check_scenario.jl", "scripts/tests/check_api_source.jl",
                      "scripts/tests/check_api_validation.jl"]),
        ("docs",     ["dev/api_plan.md", "dev/multiregion_phantom_plan.md", "docs/SCHEMA.md",
                      "docs/range_verification_recipe.md"]),
    ]
    for (group, files) in comps
        @printf("  %-9s %s\n", group * ":", join([isfile(rp(f)) ? f : "$f [MISSING]" for f in files], "\n" * " "^12))
    end
end

main()
