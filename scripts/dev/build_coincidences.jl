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

# The LOR selection core — Response, GammaAcc, pass_energy, contained_one, fill_full!,
# finish_event! — now lives in src/coincidences.jl, shared with build_true_coincidences_from_singles.jl.
# Here we only supply the CSV-row sink that `finish_event!` emits into (the truth Int8 code is
# mapped back to the "true"/"scatter" strings). The timestamps `t1,t2` [ns] and the per-pair
# residual `dt` = |Δt0| − TOF_diff are now real (an `EventTiming` is passed to finish_event!).
function write_lor_row(io, ev, x1, y1, z1, e1, t1, iz1, iphi1,
                       x2, y2, z2, e2, t2, iz2, iphi2, dt, x0, y0, z0, truth)
    ts = truth == TRUTH_TRUE ? "true" : "scatter"
    println(io, join((ev,
        round(x1, digits=4), round(y1, digits=4), round(z1, digits=4), round(e1, digits=4), round(t1, digits=4), iz1, iphi1,
        round(x2, digits=4), round(y2, digits=4), round(z2, digits=4), round(e2, digits=4), round(t2, digits=4), iz2, iphi2,
        round(dt, digits=6),
        round(x0, digits=4), round(y0, digits=4), round(z0, digits=4), ts), ","))
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(cfg_get(cfg, "output", "dir", "output")), tag)

    sigma_xyz = Float64(cfg_get(cfg, "detector", "sigma_xyz_mm", 0.0))
    eres      = Float64(cfg_get(cfg, "detector", "eres", 0.0))
    emin      = Float64(cfg_get(cfg, "detector", "emin_keV", 0.0))
    window    = Float64(cfg_get(cfg, "detector", "window_fwhm", 0.0))
    seed      = Int(cfg_get(cfg, "detector", "seed", 1234))
    fmt       = String(cfg_get(cfg, "output", "format", "csv"))
    fmt == "csv" || error("[output].format '$fmt' not supported yet (csv only; hdf5 deferred)")
    window > 0.0 && eres <= 0.0 &&
        error("[detector].window_fwhm requires eres > 0 (the window width scales with the resolution)")

    mode  = (sigma_xyz > 0.0 || eres > 0.0 || window > 0.0 || emin > 0.0) ? "det" : "truth"
    stack = isempty(a["stack"]) ? joinpath(outdir, "stack.csv") : rp(a["stack"])
    out   = isempty(a["out"])   ? joinpath(outdir, "coincidences_$mode.csv") : rp(a["out"])
    isfile(stack) || error("stack file '$stack' not found (run simulate_phantom.jl first)")

    resp = Response(sigma_xyz, eres, emin, window > 0.0, window * energy_fwhm(511.0, eres))
    rng  = MersenneTwister(seed)
    # The full stack carries no per-photon time, so (unlike the production singles) we compute it
    # here, once each gamma is complete: t = TOF(emit→first hit) + scintillation jitter (needs the
    # gamma's total deposited energy, known only at the event boundary). `stamp_t` returns it for a
    # contained gamma; the caller passes it (and the phantom-scatter flag) to finish_event!.
    crystal = String(cfg_get(cfg, "transport", "crystal_material", "CsI"))
    cryst   = load_material(joinpath(REPO, "data"), crystal)
    stamp_t(acc, ex, ey, ez) = acc.reached ?
        tof_ns((ex, ey, ez), (acc.x, acc.y, acc.z)) + first_photon_jitter(cryst, acc.e * 1e-3, rng) : 0.0

    if mode == "det"
        ecut = resp.apply_window ? "window ±$(round(resp.win_half,digits=1)) keV" :
               (emin > 0.0 ? "Emin $(round(emin,digits=0)) keV" : "no energy cut")
        println("run '$tag' [det]: σ_xyz=$(resp.sigma_xyz) mm, eres=$(round(100*resp.eres,digits=1))% @511 keV, " *
                "$ecut  (seed $seed)")
    else
        println("run '$tag' [truth]: no smearing, no energy cut")
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
        gps = [false, false]             # per-gamma phantom-scatter flag (OR-ed over the rows)
        cur_ev = -1
        ev_x0 = 0.0; ev_y0 = 0.0; ev_z0 = 0.0     # this event's emission point
        n_ev = 0; n_pair = 0; n_true = 0

        mkpath(dirname(out))
        open(out, "w") do io
            println(io, "event,x1_mm,y1_mm,z1_mm,e1_keV,t1_ns,iz1,iphi1," *
                        "x2_mm,y2_mm,z2_mm,e2_keV,t2_ns,iz2,iphi2,dt_ns,x0_mm,y0_mm,z0_mm,truth")
            for line in eachline(fin)
                isempty(line) && continue
                f = split(line, ',')
                ev = parse(Int, f[iev])
                if ev != cur_ev
                    if cur_ev != -1
                        n_ev += 1
                        t1 = stamp_t(g[1], ev_x0, ev_y0, ev_z0); t2 = stamp_t(g[2], ev_x0, ev_y0, ev_z0)
                        emitted, is_true = finish_event!((a...) -> write_lor_row(io, a...),
                                                         cur_ev, g[1], g[2], t1, t2, gps[1], gps[2],
                                                         ev_x0, ev_y0, ev_z0, resp, rng)
                        n_pair += emitted; n_true += (emitted && is_true)
                    end
                    reset!(g[1]); reset!(g[2]); gps[1] = false; gps[2] = false; cur_ev = ev
                end
                ev_x0 = parse(Float64, f[ix0]); ev_y0 = parse(Float64, f[iy0]); ev_z0 = parse(Float64, f[iz0])

                gi = parse(Int, f[ig])
                (gi == 1 || gi == 2) || continue
                acc = g[gi]
                gps[gi] |= (f[iph] == "1")
                f[ivol] == "scanner" || continue
                edep = parse(Float64, f[ide])
                edep > 0.0 || continue                  # skip the scanner :escape row (e_dep=0)
                fill_full!(acc, parse(Float64, f[ix]), parse(Float64, f[iy]), parse(Float64, f[iz_]),
                           edep, parse(Int, f[iiz]), parse(Int, f[iiphi]))
            end
            if cur_ev != -1                              # the final event
                n_ev += 1
                t1 = stamp_t(g[1], ev_x0, ev_y0, ev_z0); t2 = stamp_t(g[2], ev_x0, ev_y0, ev_z0)
                emitted, is_true = finish_event!((a...) -> write_lor_row(io, a...),
                                                 cur_ev, g[1], g[2], t1, t2, gps[1], gps[2],
                                                 ev_x0, ev_y0, ev_z0, resp, rng)
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
