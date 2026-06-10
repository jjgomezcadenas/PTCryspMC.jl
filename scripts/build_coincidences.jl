#!/usr/bin/env julia
# Build a TRUTH-level list-mode coincidence file from a navigated back-to-back stack
# (scripts/navigate_back_to_back.jl output). One record per accepted gamma pair.
#
# Selection (truth, no detector resolution yet):
#   • both gammas of the annihilation interact in the detector (≥1 scanner deposit each), and
#   • each gamma is contained in exactly ONE crystal block (no overspill to a neighbour;
#     partial energy is allowed — a gamma that escapes the detector from a single crystal
#     still gives a clean line of response).
#
# Per gamma the record carries the first-interaction point (the LOR point), the summed
# energy in that crystal, the block (iz, iphi), and a time t — a DUMMY 0.0 for now: the
# simulation has no timing yet, and real times are assigned by the randoms pass (Step 5).
# The pair is tagged `true` if neither gamma Compton-scattered in the phantom, else `scatter`.
#
# This is the same-annihilation, truth-level skeleton of coincidences_<config>.csv; the
# detector smearing (σ_xyz, σ_t, FWHM(E)) and the ±2·FWHM energy window are Step 4.
#
# Streaming: the stack is event-ordered, so we read line by line, accumulate one event at
# a time, emit its coincidence when the event number changes, and discard — O(1) memory,
# no whole-file load (the point of doing this in Julia for large stacks).
#
# Run from the repo root:
#   julia --project=. scripts/build_coincidences.jl --stack output/nav_b2b_bgo_stack.csv

using ArgParse

function parse_cli()
    s = ArgParseSettings(description="Build a truth list-mode coincidence file from a navigated back-to-back stack.")
    @add_arg_table! s begin
        "--stack"; help = "input stack CSV (navigate_back_to_back.jl output)"; required = true
        "--out";   help = "output coincidence CSV (default output/coincidences_<stack>.csv)"; default = ""
    end
    parse_args(s)
end

# Per-gamma accumulator, reset at each event boundary. We never store the row list: the
# first scanner deposit (rows are step-ordered) fixes the LOR point and the block; later
# deposits add energy and, if in a different block, flag overspill.
mutable struct GammaAcc
    reached::Bool
    x::Float64; y::Float64; z::Float64
    iz::Int; iphi::Int
    e::Float64
    overspill::Bool          # a deposit landed in a second, different block
    phscat::Bool
end
GammaAcc() = GammaAcc(false, 0.0, 0.0, 0.0, -1, -1, 0.0, false, false)

function reset!(a::GammaAcc)
    a.reached = false; a.x = a.y = a.z = 0.0; a.iz = a.iphi = -1
    a.e = 0.0; a.overspill = false; a.phscat = false
    a
end

"A gamma passes if it reached the detector and stayed within one crystal block."
contained_one(a::GammaAcc) = a.reached && !a.overspill

# Finalise one event's two gammas: emit a coincidence row if both pass. Returns
# (emitted::Bool, is_true::Bool).
function finish_event!(io, ev::Int, g1::GammaAcc, g2::GammaAcc)
    (contained_one(g1) && contained_one(g2)) || return (false, false)
    is_true = !(g1.phscat || g2.phscat)
    truth = is_true ? "true" : "scatter"
    println(io, join((ev,
        round(g1.x, digits=4), round(g1.y, digits=4), round(g1.z, digits=4),
        round(g1.e, digits=4), 0.0, g1.iz, g1.iphi,
        round(g2.x, digits=4), round(g2.y, digits=4), round(g2.z, digits=4),
        round(g2.e, digits=4), 0.0, g2.iz, g2.iphi, truth), ","))
    (true, is_true)
end

function main()
    a = parse_cli()
    stack = a["stack"]
    isfile(stack) || error("stack file '$stack' not found")
    out = isempty(a["out"]) ?
        joinpath(dirname(stack), "coincidences_" * basename(stack)) : a["out"]

    open(stack, "r") do fin
        header = split(strip(readline(fin)), ',')
        col = Dict(String(h) => i for (i, h) in enumerate(header))
        for c in ("event_number", "gamma", "x_mm", "y_mm", "z_mm", "e_dep_keV",
                  "volume", "iz", "iphi", "phantom_scatter")
            haskey(col, c) || error("stack is missing column '$c'")
        end
        iev = col["event_number"]; ig = col["gamma"]
        ix = col["x_mm"]; iy = col["y_mm"]; iz_ = col["z_mm"]; ide = col["e_dep_keV"]
        ivol = col["volume"]; iiz = col["iz"]; iiphi = col["iphi"]; iph = col["phantom_scatter"]

        g = (GammaAcc(), GammaAcc())     # gamma 1, gamma 2
        cur_ev = -1
        n_ev = 0; n_pair = 0; n_true = 0

        mkpath(dirname(out))
        open(out, "w") do io
            println(io, "event,x1_mm,y1_mm,z1_mm,e1_keV,t1_ns,iz1,iphi1," *
                        "x2_mm,y2_mm,z2_mm,e2_keV,t2_ns,iz2,iphi2,truth")
            for line in eachline(fin)
                isempty(line) && continue
                f = split(line, ',')
                ev = parse(Int, f[iev])
                if ev != cur_ev
                    if cur_ev != -1
                        n_ev += 1
                        emitted, is_true = finish_event!(io, cur_ev, g[1], g[2])
                        n_pair += emitted; n_true += (emitted && is_true)
                    end
                    reset!(g[1]); reset!(g[2]); cur_ev = ev
                end

                gi = parse(Int, f[ig])
                (gi == 1 || gi == 2) || continue
                acc = g[gi]
                acc.phscat |= (f[iph] == "1")
                f[ivol] == "scanner" || continue
                edep = parse(Float64, f[ide])
                edep > 0.0 || continue                  # skip the scanner :escape row (e_dep=0)
                blk = (parse(Int, f[iiz]), parse(Int, f[iiphi]))
                if !acc.reached                          # first scanner deposit: the LOR point
                    acc.reached = true
                    acc.x = parse(Float64, f[ix]); acc.y = parse(Float64, f[iy]); acc.z = parse(Float64, f[iz_])
                    acc.iz, acc.iphi = blk
                elseif blk != (acc.iz, acc.iphi)
                    acc.overspill = true                 # a second, different crystal
                end
                acc.e += edep
            end
            if cur_ev != -1                              # the final event
                n_ev += 1
                emitted, is_true = finish_event!(io, cur_ev, g[1], g[2])
                n_pair += emitted; n_true += (emitted && is_true)
            end
        end

        println("read $n_ev events from $stack")
        println("wrote $n_pair coincidences -> $out")
        if n_pair > 0
            pct = round(100 * n_pair / max(n_ev, 1), digits=1)
            tpct = round(100 * n_true / n_pair, digits=1)
            println("  coincidence acceptance: $pct% of annihilations")
            println("  truth split: $n_true true ($tpct%), $(n_pair - n_true) scatter ($(round(100 - tpct, digits=1))%)")
        end
    end
end

main()
