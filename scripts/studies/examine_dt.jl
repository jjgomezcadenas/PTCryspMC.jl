#!/usr/bin/env julia
# Examine the timing of the TRUE coincidences to choose the coincidence window τ.
#
# Streams a truth LOR file (build_true_coincidences_from_singles.jl output) and dumps, per LOR,
# the two timing quantities:
#   • dT_ns = |t1 − t2|         — the RAW timestamp difference; THIS is the coincidence-window
#                                 variable (a scanner windows on raw arrival-time difference, and
#                                 it is the only quantity also defined for randoms).
#   • dt_ns = |Δt0| − TOF_diff  — the TOF-corrected residual = the timing RESOLUTION (diagnostic;
#                                 undefined for randoms, so NOT the cut).
# Plot dT_ns (py/plot_dt.py) and read the window τ off its tail; the reco stage then keeps trues
# and randoms with dT ≤ τ. Every truth LOR is a genuine same-annihilation coincidence at this
# stage (true + scatter share the same timing), so no truth split is applied.
#
#   julia --project=. scripts/studies/examine_dt.jl --config runs/sphere_water_csi.toml
#   julia --project=. scripts/studies/examine_dt.jl --config runs/sphere_water_csi.toml --lors path.h5

using PTCryspMC
using ArgParse
using Printf

function parse_cli()
    s = ArgParseSettings(description="Dump dT/dt of the true coincidences for the τ decision.")
    @add_arg_table! s begin
        "--config"; help = "run config TOML"; required = true
        "--lors";   help = "override the LOR path (default prod/<tag>/lors_truth.h5)"; default = ""
        "--out";    help = "override the dt.csv path (default studies/dt/<tag>_dt.csv)"; default = ""
    end
    parse_args(s)
end

# Quantiles from an (unsorted) vector, sorting once. q ∈ [0,1].
function quantiles(v::Vector{Float32}, qs)
    sort!(v); n = length(v)
    [n == 0 ? NaN32 : v[clamp(ceil(Int, q * n), 1, n)] for q in qs]
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    tag = run_tag(cfg, a["config"])
    lors = isempty(a["lors"]) ? joinpath(rp(prod_base(cfg)), tag, "lors_truth.h5") : rp(a["lors"])
    isfile(lors) || error("LOR file '$lors' not found (run build_true_coincidences_from_singles.jl first)")
    out = isempty(a["out"]) ? joinpath(REPO, "studies", "dt", "$(tag)_dt.csv") : rp(a["out"])
    mkpath(dirname(out))

    println("examine dt: $lors")
    dTs = Float32[]                      # collected for the summary (CSV is written streaming)
    n = 0
    open(out, "w") do io
        println(io, "dT_ns,dt_ns")
        foreach_coincidences_hdf5(lors) do b
            for i in 1:length(b)
                dT = abs(b.t1[i] - b.t2[i])
                println(io, @sprintf("%.5f,%.6f", dT, b.dt[i]))
                push!(dTs, dT); n += 1
            end
        end
    end

    p50, p99, p999 = quantiles(dTs, (0.5, 0.99, 0.999))
    mx = isempty(dTs) ? NaN32 : dTs[end]                 # sorted by quantiles()
    mn = isempty(dTs) ? NaN32 : sum(dTs) / n
    println("wrote $n rows -> $out")
    @printf("  dT_ns = |t1-t2| (the coincidence-window variable):\n")
    @printf("    mean %.4f   p50 %.4f   p99 %.4f   p99.9 %.4f   max %.4f  [ns]\n", mn, p50, p99, p999, mx)
    println("  → plot py/plot_dt.py and read τ off the dT tail (trade trues kept vs randoms ∝ τ).")
end

main()
