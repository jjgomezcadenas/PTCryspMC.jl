#!/usr/bin/env julia
# Build the TRUE list-mode LOR file from a SINGLES stack (scripts/simulate_source_mt.jl output) —
# the same-annihilation coincidences (true + scatter), with NO detector smearing and NO energy
# cut: ground-truth positions / energies, plus the per-gamma timestamps `t1,t2` and the residual
# `dt = |Δt0| − TOF_diff`. One record per accepted coincidence = one LOR. The singles row already
# IS the formed per-gamma hit, so this fills `GammaAcc` directly and runs the shared selection
# core (src/coincidences.jl).
#
# This is the clean input to the DT study: examine the `t1,t2,dt` distribution (examine_dt.jl) to
# choose the coincidence window τ. Detector smearing, the DT cut and the randoms come AFTER τ is
# chosen, in the reco stage → lors_det.h5 (not produced here).
#
# Input: singles in either format (default prod/<tag>/singles.{h5,csv} by [output].format;
# --singles override; format by extension). Output: prod/<tag>/lors_truth.h5 — quantized Int16
# columns, chunked + shuffle+deflate, truth ∈ {true=0, scatter=1}, `has_randoms=false`.
#
# Run from the repo root:
#   julia --project=. scripts/build_true_coincidences_from_singles.jl --config runs/sphere_water_csi.toml

using PTCryspMC
using ArgParse
using Random

function parse_cli()
    s = ArgParseSettings(description="Build the LOR list from a singles stack (TOML-config driven).")
    @add_arg_table! s begin
        "--config";  help = "run config TOML"; required = true
        "--singles"; help = "override the input singles path (default prod/<tag>/singles.{h5,csv})"; default = ""
        "--out";     help = "override the output LOR path"; default = ""
    end
    parse_args(s)
end

# Streaming event-grouping state (persists across HDF5 batches).
mutable struct LorState
    cur_ev::Int; maxev::Int
    ev_x0::Float64; ev_y0::Float64; ev_z0::Float64
    ev_tdec::Float64                 # this event's absolute decay time [s] (event_time)
    n_pair::Int; n_true::Int
end
LorState() = LorState(-1, 0, 0.0, 0.0, 0.0, 0.0, 0, 0)

# Feed one singles row: close the previous event at a boundary (finish_event! → writer), then
# fill the gamma's accumulator directly. `g` = (GammaAcc, GammaAcc).
# `g` = (GammaAcc, GammaAcc); `gt`/`gps` carry the per-gamma timestamp / phantom-scatter count (not
# stored on the accumulator) until finish_event! at the event boundary. `acts` is the per-isotope
# activity-model vector (API) or length-1 (clinic/count, isotope always 0); the row's isotope picks
# the model whose event_time stamps the event's t_decay (finish_event! is decay-time-agnostic —
# the sink closure appends it).
function feed_row!(st::LorState, g, gt, gps, w, resp, rng, acts, ev::Int, gi::Int, iso::Int,
                   x::Float64, y::Float64, z::Float64, e::Float64, iz::Int, iphi::Int,
                   nblocks::Int, nscat::Int, t_rel::Float64, x0::Float64, y0::Float64, z0::Float64)
    if ev != st.cur_ev
        if st.cur_ev != -1
            emitted, is_true = finish_event!((a...) -> push_coincidence!(w, a..., st.ev_tdec),
                                             st.cur_ev, g[1], g[2], gt[1], gt[2], gps[1], gps[2],
                                             st.ev_x0, st.ev_y0, st.ev_z0, resp, rng)
            st.n_pair += emitted; st.n_true += (emitted && is_true)
        end
        reset!(g[1]); reset!(g[2]); gt[1] = 0.0; gt[2] = 0.0; gps[1] = 0; gps[2] = 0; st.cur_ev = ev
    end
    ev > st.maxev && (st.maxev = ev)
    st.ev_x0 = x0; st.ev_y0 = y0; st.ev_z0 = z0
    st.ev_tdec = event_time(acts[iso + 1], ev)
    (gi == 1 || gi == 2) || return
    fill_singles!(g[gi], x, y, z, e, iz, iphi, nblocks)
    gt[gi] = t_rel; gps[gi] = nscat
    return
end

function feed_singles_csv!(st, g, gt, gps, w, resp, rng, acts, path)
    open(path, "r") do io
        header = split(strip(readline(io)), ',')
        col = Dict(String(h) => i for (i, h) in enumerate(header))
        for c in ("event_number", "gamma", "x_mm", "y_mm", "z_mm", "e_keV", "iz", "iphi",
                  "nblocks", "n_scatter", "isotope", "x0_mm", "y0_mm", "z0_mm", "t_rel_ns")
            haskey(col, c) || error("singles stack is missing column '$c'")
        end
        ie = col["event_number"]; ig = col["gamma"]; ix = col["x_mm"]; iy = col["y_mm"]; iz = col["z_mm"]
        iee = col["e_keV"]; iiz = col["iz"]; iip = col["iphi"]; inb = col["nblocks"]; iph = col["n_scatter"]
        iis = col["isotope"]; ix0 = col["x0_mm"]; iy0 = col["y0_mm"]; iz0 = col["z0_mm"]; it = col["t_rel_ns"]
        for line in eachline(io)
            isempty(line) && continue
            f = split(line, ',')
            feed_row!(st, g, gt, gps, w, resp, rng, acts, parse(Int, f[ie]), parse(Int, f[ig]),
                      parse(Int, f[iis]),
                      parse(Float64, f[ix]), parse(Float64, f[iy]), parse(Float64, f[iz]),
                      parse(Float64, f[iee]), parse(Int, f[iiz]), parse(Int, f[iip]),
                      parse(Int, f[inb]), parse(Int, f[iph]), parse(Float64, f[it]),
                      parse(Float64, f[ix0]), parse(Float64, f[iy0]), parse(Float64, f[iz0]))
        end
    end
end

function feed_singles_hdf5!(st, g, gt, gps, w, resp, rng, acts, path)
    foreach_singles_hdf5(path) do b
        for i in 1:length(b)
            feed_row!(st, g, gt, gps, w, resp, rng, acts, Int(b.event[i]), Int(b.gamma[i]),
                      Int(b.isotope[i]),
                      decode_xyz(b.x[i]), decode_xyz(b.y[i]), decode_xyz(b.z[i]), decode_e(b.e[i]),
                      Int(b.iz[i]), Int(b.iphi[i]), Int(b.nblocks[i]), Int(b.n_scatter[i]),
                      Float64(b.t_rel[i]),
                      decode_xyz(b.x0[i]), decode_xyz(b.y0[i]), decode_xyz(b.z0[i]))
        end
    end
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(prod_base(cfg)), tag)

    # Per-photon timestamps `t_rel` are already stamped in the singles; this just reads them. The
    # activity model stamps each event's absolute decay time `t_decay_s` (event_time(ev), the same
    # per-isotope clock the randoms pass uses to restore the absolute time). Crystal name is
    # recorded for the attrs.
    crystal = String(cfg_get(cfg, "transport", "crystal_material", "CsI"))
    act     = ActivityModel(cfg)

    # Truth coincidences: no smearing, no energy cut (Response all-off; rng unused but kept for the
    # shared finish_event! signature).
    seed = Int(cfg_get(cfg, "transport", "seed", 1234))
    resp = Response(0.0, 0.0, 0.0, false, 0.0)

    fmt = String(cfg_get(cfg, "output", "format", "csv"))
    default_in = joinpath(outdir, fmt == "hdf5" ? "singles.h5" : "singles.csv")
    singles = isempty(a["singles"]) ? default_in : rp(a["singles"])
    isfile(singles) || error("singles file '$singles' not found (run simulate_source_mt.jl first)")
    ishdf5 = endswith(singles, ".h5")
    out = isempty(a["out"]) ? joinpath(outdir, "lors_truth.h5") : rp(a["out"])

    # Decay-time clock for t_decay_s: per-isotope when the singles carry `isotope_half_lives`
    # (API — one truncated-exponential model per isotope over [0, t_meas]), else the single config
    # `ActivityModel` (clinic/count; the isotope column is all 0, so acts[1] is used for every
    # single). Identical construction to build_randoms_from_singles — the two passes MUST share
    # the clock so t_decay_s and the randoms pairing agree.
    hls = ishdf5 ? Vector{Float64}(singles_hdf5_attr(singles, "isotope_half_lives", Float64[])) : Float64[]
    if !isempty(hls)
        tmeas = Float64(singles_hdf5_attr(singles, "t_meas_s", 600.0))
        tseed = Int(singles_hdf5_attr(singles, "time_seed", 1234))
        acts  = ActivityModel[ActivityModel(; t0=0.0, t1=tmeas, half_life_s=h, seed=tseed) for h in hls]
    else
        acts  = ActivityModel[act]
    end

    println("run '$tag' [truth]: same-annihilation coincidences, no smearing/cut")
    println("  crystal: $crystal  |  singles ($(ishdf5 ? "hdf5" : "csv")): $singles")

    meta = Dict{String,Any}("scenario_tag" => tag, "mode" => "truth", "has_randoms" => false,
        "crystal" => crystal, "seed" => seed, "t_relative_to_decay" => true,
        "t0_s" => acts[1].t0, "t1_s" => acts[1].t1, "n_isotopes" => length(acts),
        "half_life_s" => log(2.0) / act.λ, "time_seed" => Int(acts[1].seed))
    w = CoincidenceWriter(out, meta)
    st = LorState()
    g = (GammaAcc(), GammaAcc())
    gt = [0.0, 0.0]; gps = [0, 0]                    # per-gamma timestamp / phantom-scatter count
    rng = MersenneTwister(seed)

    ishdf5 ? feed_singles_hdf5!(st, g, gt, gps, w, resp, rng, acts, singles) :
             feed_singles_csv!(st, g, gt, gps, w, resp, rng, acts, singles)
    if st.cur_ev != -1                                   # the final event
        emitted, is_true = finish_event!((a...) -> push_coincidence!(w, a..., st.ev_tdec),
                                         st.cur_ev, g[1], g[2], gt[1], gt[2], gps[1], gps[2],
                                         st.ev_x0, st.ev_y0, st.ev_z0, resp, rng)
        st.n_pair += emitted; st.n_true += (emitted && is_true)
    end
    # Total annihilations for the acceptance denominator (the singles stack has no rows for
    # both-missed events). HDF5 carries the exact `nevents` attribute; CSV has no metadata, so
    # use the max event seen (≈ nevents). Stamp it onto the LOR file.
    nevents = ishdf5 ? Int(singles_hdf5_attr(singles, "nevents", 0)) : 0
    nevents > 0 || (nevents = st.maxev)
    set_lor_attr!(w, "nevents", nevents)
    # Transport provenance for exact singles regeneration after pruning: the actual transport seed
    # and chunk count (with nevents + the config, the full recipe). Carried from the singles attrs;
    # CSV singles have no metadata, so they stay absent there.
    if ishdf5
        set_lor_attr!(w, "transport_seed", Int(singles_hdf5_attr(singles, "seed", seed)))
        set_lor_attr!(w, "nchunks", Int(singles_hdf5_attr(singles, "nchunks", 0)))
        copy_provenance!(w, singles)     # carry the source/geometry provenance (scenario, dose, …)
    end
    nlor = close(w)

    println("wrote $nlor LORs -> $out")
    if nlor > 0
        denom = nevents > 0 ? nevents : st.n_pair
        tpct = round(100 * st.n_true / nlor, digits=1)
        println("  acceptance: $(round(100*nlor/max(denom,1), digits=1))% of $nevents annihilations")
        println("  truth split: $(st.n_true) true ($tpct%), $(nlor - st.n_true) scatter ($(round(100-tpct, digits=1))%)")
        println("  (truth ∈ {true,scatter}; randoms added later by reco_lors; has_randoms=false)")
    end
end

main()
