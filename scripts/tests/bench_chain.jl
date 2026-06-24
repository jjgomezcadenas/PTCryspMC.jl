#!/usr/bin/env julia
# Benchmark the production LOR chain on an existing singles stack — wall time + peak RSS per stage,
# to decide where MT / streaming is needed for the 1e8 production. Each stage runs as a subprocess
# under `/usr/bin/time -l` (macOS); we parse its "maximum resident set size". Subtracting a measured
# bare-julia-load baseline gives the compute time, which we extrapolate to `--target` events
# (streaming stages ~linear in N; build_randoms holds an in-memory store ~linear in N → the RAM to
# watch). Run after the singles exist (simulate_source_mt.jl), e.g.:
#
#   julia --project=. scripts/tests/bench_chain.jl --config runs/sphere_water_csi.toml \
#       --singles prod/sphere_water_csi/singles.h5 --target 100000000

using PTCryspMC
using ArgParse
using Printf

function parse_cli()
    s = ArgParseSettings(description="Benchmark the LOR chain (wall + peak RSS per stage).")
    @add_arg_table! s begin
        "--config";  help = "run config TOML"; required = true
        "--singles"; help = "singles stack (default prod/<tag>/singles.h5)"; default = ""
        "--target";  help = "events to extrapolate to"; arg_type = Int; default = 100_000_000
    end
    parse_args(s)
end

# Run `cmd` under /usr/bin/time -l; return (wall_s, peak_rss_bytes). Stage stdout is muted; its
# stderr + the time report are captured (the RSS line lives there).
function timed(cmd::Cmd)
    err = IOBuffer()
    t0 = time()
    ok = try
        run(pipeline(`/usr/bin/time -l $cmd`; stdout = devnull, stderr = err)); true
    catch
        false
    end
    wall = time() - t0
    txt = String(take!(err))
    m = match(r"(\d+)\s+maximum resident set size", txt)
    (wall, m === nothing ? 0 : parse(Int, m.captures[1]), ok)
end

gb(bytes) = bytes / 2^30

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    cfg = read_config(a["config"]); tag = run_tag(cfg, a["config"])
    singles = isempty(a["singles"]) ? joinpath(REPO, prod_base(cfg), tag, "singles.h5") : a["singles"]
    isfile(singles) || error("singles '$singles' not found (run simulate_source_mt.jl first)")
    N = Int(singles_hdf5_attr(singles, "nevents", 0))
    N > 0 || error("singles file has no nevents attribute")
    target = a["target"]

    jl  = `julia --project=$REPO`
    scr(p) = joinpath(REPO, "scripts", p)
    # :flat  = streaming, RAM ~constant in N (HDF5 batch + fixed buffers);
    # :linear = holds an in-memory store ~N (the τ-bin / external-sort candidate).
    stages = [
        ("build_true",    `$jl $(scr("build_true_coincidences_from_singles.jl")) --config $(a["config"]) --singles $singles`, :flat),
        ("build_randoms", `$jl $(scr("build_randoms_from_singles.jl")) --config $(a["config"]) --singles $singles`, :linear),
        ("reco",          `$jl $(scr("reco_lors.jl")) --config $(a["config"])`, :flat),
    ]

    println("benchmark '$tag' on $N events  ($singles)")
    print("  measuring bare-julia load baseline… "); flush(stdout)
    base_wall, base_rss, _ = timed(`$jl -e "using PTCryspMC"`)
    @printf "%.2fs, %.2f GB\n\n" base_wall gb(base_rss)

    @printf "  %-14s %6s %8s %9s %9s | %10s %10s\n" "stage" "RAM" "wall[s]" "comp[s]" "RSS[GB]" "→$(target)[s]" "→RAM[GB]"
    results = NamedTuple[]
    for (name, cmd, scaling) in stages
        wall, rss, ok = timed(cmd)
        ok || @printf "  %-14s  FAILED\n" name
        ok || continue
        compute = max(wall - base_wall, 1e-3)
        scale = target / N
        t_target = compute * scale                          # time ~linear in N for every stage
        # RAM: :flat stays put (streaming); :linear grows as base + data·N.
        ram_target = scaling === :linear ? gb(base_rss) + gb(max(rss - base_rss, 0)) * scale : gb(rss)
        @printf "  %-14s %6s %8.2f %9.2f %9.2f | %10.0f %10.1f\n" name String(scaling) wall compute gb(rss) t_target ram_target
        push!(results, (; name, compute, rss, t_target, ram_target))
    end

    println()
    tot = sum(r.t_target for r in results; init = 0.0)
    @printf "  extrapolated chain compute @ %d events: %.0f s (%.1f min), serial\n" target tot tot/60
    slow = isempty(results) ? nothing : argmax(r -> r.t_target, results)
    big  = isempty(results) ? nothing : argmax(r -> r.ram_target, results)
    slow === nothing || @printf "  slowest stage: %s (%.0f s) — MT candidate if it dominates\n" slow.name slow.t_target
    big  === nothing || @printf "  largest RAM:   %s (%.1f GB) — the streaming/τ-bin candidate\n" big.name big.ram_target
    println("  (rough: wall includes ~$(round(base_wall,digits=1))s julia load; time ~linear, RAM = base + data·N)")
end

main()
