#!/usr/bin/env julia
# Build a list-mode coincidence file from a phantom-simulation stack
# (scripts/simulate_phantom.jl output). One record per accepted gamma pair.
#
# Selection: a pair is accepted when
#   • both gammas of the annihilation interact in the detector (≥1 scanner deposit each), and
#   • each gamma is contained in exactly ONE crystal block (no overspill to a neighbour;
#     partial energy is allowed — a gamma that escapes the detector from a single crystal
#     still gives a clean line of response), and
#   • (detector mode only) both smeared energies fall inside the energy window.
#
# Detector response comes from the run config's [detector] section (all OFF → TRUTH mode):
#   sigma_xyz_mm  Gaussian smear of each hit position [mm] (includes DOI)
#   eres          energy resolution: fractional FWHM at 511 keV; FWHM(E)=a·√(511/E)
#   window_fwhm   energy window half-width in FWHM(511) units → keep |E−511| ≤ w·a·511 keV
#   seed          RNG seed for the smearing
#
# Per gamma the record carries the first-interaction point (the LOR point, smeared in
# detector mode), the summed energy in that crystal (smeared), the block (iz,iphi), and a
# time t — a DUMMY 0.0 (the simulation has no timing yet; real times come from the randoms
# pass, Step 5). The true emission point (x0,y0,z0) is carried through for source
# validation. The pair is tagged `true` if neither gamma scattered in the phantom, else
# `scatter` (smearing does not change the truth tag).
#
# TOML-config driven: the stack is read from output/<tag>/stack.csv and the coincidences
# written to output/<tag>/coincidences_{truth|det}.csv (tag = config filename). Streaming:
# the stack is event-ordered, so we read line by line, accumulate one event at a time, emit
# its coincidence when the event number changes, and discard — O(1) memory.
#
# Run from the repo root:
#   julia --project=. scripts/build_coincidences.jl --config runs/sphere_water_csi.toml

using PTCryspMC                       # detector response + config: smear_*, energy_fwhm, read_config
using ArgParse
using Random

function parse_cli()
    s = ArgParseSettings(description="Build a list-mode coincidence file from a phantom-simulation stack (TOML-config driven).")
    @add_arg_table! s begin
        "--config"; help = "run config TOML"; required = true
        "--stack";  help = "override the input stack path (default output/<tag>/stack.csv)"; default = ""
        "--out";    help = "override the output coincidence path"; default = ""
    end
    parse_args(s)
end

# Detector response. All-zero / no-window = truth mode (no smearing, no cut).
struct Response
    sigma_xyz::Float64       # mm
    eres::Float64            # fractional FWHM at 511 keV
    apply_window::Bool
    win_half::Float64        # window half-width [keV]
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

# Finalise one event's two gammas: smear, apply the energy window, and emit a coincidence
# row if both pass. (x0,y0,z0) is the true emission point. Returns (emitted, is_true).
function finish_event!(io, ev::Int, g1::GammaAcc, g2::GammaAcc,
                       x0::Float64, y0::Float64, z0::Float64, resp::Response, rng)
    (contained_one(g1) && contained_one(g2)) || return (false, false)
    e1 = smear_energy(g1.e, resp.eres, rng)
    e2 = smear_energy(g2.e, resp.eres, rng)
    if resp.apply_window
        (abs(e1 - 511.0) <= resp.win_half && abs(e2 - 511.0) <= resp.win_half) || return (false, false)
    end
    x1, y1, z1 = smear_position((g1.x, g1.y, g1.z), resp.sigma_xyz, rng)
    x2, y2, z2 = smear_position((g2.x, g2.y, g2.z), resp.sigma_xyz, rng)
    is_true = !(g1.phscat || g2.phscat)
    truth = is_true ? "true" : "scatter"
    println(io, join((ev,
        round(x1, digits=4), round(y1, digits=4), round(z1, digits=4), round(e1, digits=4), 0.0, g1.iz, g1.iphi,
        round(x2, digits=4), round(y2, digits=4), round(z2, digits=4), round(e2, digits=4), 0.0, g2.iz, g2.iphi,
        round(x0, digits=4), round(y0, digits=4), round(z0, digits=4), truth), ","))
    (true, is_true)
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(cfg_get(cfg, "output", "dir", "output")), tag)

    sigma_xyz = Float64(cfg_get(cfg, "detector", "sigma_xyz_mm", 0.0))
    eres      = Float64(cfg_get(cfg, "detector", "eres", 0.0))
    window    = Float64(cfg_get(cfg, "detector", "window_fwhm", 0.0))
    seed      = Int(cfg_get(cfg, "detector", "seed", 1234))
    fmt       = String(cfg_get(cfg, "output", "format", "csv"))
    fmt == "csv" || error("[output].format '$fmt' not supported yet (csv only; hdf5 deferred)")
    window > 0.0 && eres <= 0.0 &&
        error("[detector].window_fwhm requires eres > 0 (the window width scales with the resolution)")

    mode  = (sigma_xyz > 0.0 || eres > 0.0 || window > 0.0) ? "det" : "truth"
    stack = isempty(a["stack"]) ? joinpath(outdir, "stack.csv") : rp(a["stack"])
    out   = isempty(a["out"])   ? joinpath(outdir, "coincidences_$mode.csv") : rp(a["out"])
    isfile(stack) || error("stack file '$stack' not found (run simulate_phantom.jl first)")

    resp = Response(sigma_xyz, eres, window > 0.0, window * energy_fwhm(511.0, eres))
    rng  = MersenneTwister(seed)

    if mode == "det"
        println("run '$tag' [det]: σ_xyz=$(resp.sigma_xyz) mm, eres=$(round(100*resp.eres,digits=1))% @511 keV, " *
                (resp.apply_window ? "window ±$(round(resp.win_half,digits=1)) keV" : "no window") *
                "  (seed $seed)")
    else
        println("run '$tag' [truth]: no smearing, no energy window")
    end

    open(stack, "r") do fin
        header = split(strip(readline(fin)), ',')
        col = Dict(String(h) => i for (i, h) in enumerate(header))
        for c in ("event_number", "gamma", "x_mm", "y_mm", "z_mm", "e_dep_keV",
                  "volume", "iz", "iphi", "phantom_scatter", "x0_mm", "y0_mm", "z0_mm")
            haskey(col, c) || error("stack is missing column '$c'")
        end
        iev = col["event_number"]; ig = col["gamma"]
        ix = col["x_mm"]; iy = col["y_mm"]; iz_ = col["z_mm"]; ide = col["e_dep_keV"]
        ivol = col["volume"]; iiz = col["iz"]; iiphi = col["iphi"]; iph = col["phantom_scatter"]
        ix0 = col["x0_mm"]; iy0 = col["y0_mm"]; iz0 = col["z0_mm"]

        g = (GammaAcc(), GammaAcc())     # gamma 1, gamma 2
        cur_ev = -1
        ev_x0 = 0.0; ev_y0 = 0.0; ev_z0 = 0.0     # this event's emission point
        n_ev = 0; n_pair = 0; n_true = 0

        mkpath(dirname(out))
        open(out, "w") do io
            println(io, "event,x1_mm,y1_mm,z1_mm,e1_keV,t1_ns,iz1,iphi1," *
                        "x2_mm,y2_mm,z2_mm,e2_keV,t2_ns,iz2,iphi2,x0_mm,y0_mm,z0_mm,truth")
            for line in eachline(fin)
                isempty(line) && continue
                f = split(line, ',')
                ev = parse(Int, f[iev])
                if ev != cur_ev
                    if cur_ev != -1
                        n_ev += 1
                        emitted, is_true = finish_event!(io, cur_ev, g[1], g[2], ev_x0, ev_y0, ev_z0, resp, rng)
                        n_pair += emitted; n_true += (emitted && is_true)
                    end
                    reset!(g[1]); reset!(g[2]); cur_ev = ev
                end
                ev_x0 = parse(Float64, f[ix0]); ev_y0 = parse(Float64, f[iy0]); ev_z0 = parse(Float64, f[iz0])

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
                emitted, is_true = finish_event!(io, cur_ev, g[1], g[2], ev_x0, ev_y0, ev_z0, resp, rng)
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
