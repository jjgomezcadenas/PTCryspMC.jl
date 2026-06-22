#!/usr/bin/env julia
# Example / illustration of the activity model (src/activity.jl): sample the annihilation times
# of many events from the ¹⁵O decay over an acquisition window and histogram them, so the
# front-loaded truncated-exponential distribution can be plotted against the analytic curve
# (py/plot_o15_lifetime.py). This is the per-event time the randoms pass will draw.
#
#   julia --project=. scripts/studies/o15_lifetime.jl                       # ¹⁵O, [0, 600] s
#   julia --project=. scripts/studies/o15_lifetime.jl --t1 1200 --nevents 5000000

using PTCryspMC
using ArgParse

function parse_cli()
    s = ArgParseSettings(description="Sample + histogram the activity-model event times (example).")
    @add_arg_table! s begin
        "--half-life-s"; help = "isotope half-life [s]"; arg_type = Float64; default = O15_HALFLIFE_S
        "--t0";          help = "acquisition start [s]"; arg_type = Float64; default = 0.0
        "--t1";          help = "acquisition end [s]";   arg_type = Float64; default = 600.0
        "--nevents";     help = "events to sample";      arg_type = Int;     default = 1_000_000
        "--nbins";       help = "histogram bins";        arg_type = Int;     default = 60
        "--seed";        help = "activity time seed";    arg_type = Int;     default = 1234
        "--out";         help = "output CSV";            default = "studies/lifetime/o15_lifetime.csv"
    end
    parse_args(s)
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    out = isabspath(a["out"]) ? a["out"] : joinpath(REPO, a["out"])
    mkpath(dirname(out))

    m = ActivityModel(; t0=a["t0"], t1=a["t1"], half_life_s=a["half-life-s"], seed=a["seed"])
    N = a["nevents"]; nb = a["nbins"]
    t0 = m.t0; t1 = m.t1; w = (t1 - t0) / nb

    counts = zeros(Int, nb)
    sumt = 0.0
    for ev in 1:N
        t = event_time(m, ev)
        sumt += t
        counts[clamp(floor(Int, (t - t0) / w) + 1, 1, nb)] += 1
    end

    open(out, "w") do io
        println(io, "# o15_lifetime: half_life_s=$(a["half-life-s"]) t0_s=$t0 t1_s=$t1 " *
                    "nevents=$N lambda_per_s=$(m.λ) seed=$(a["seed"])")
        println(io, "t_lo_s,t_hi_s,count")
        for b in 1:nb
            println(io, join((round(t0 + (b-1)*w, digits=4), round(t0 + b*w, digits=4), counts[b]), ","))
        end
    end

    mean_analytic = t0 + 1/m.λ - (t1-t0)*exp(-m.λ*(t1-t0)) / (1 - exp(-m.λ*(t1-t0)))
    println("¹⁵O-style: T½=$(round(a["half-life-s"],digits=1)) s, window [$t0, $t1] s, $N events")
    println("  mean event time: sampled $(round(sumt/N, digits=1)) s vs analytic $(round(mean_analytic, digits=1)) s")
    println("  wrote $nb-bin histogram -> $out   (plot: py/plot_o15_lifetime.py)")
end

main()
