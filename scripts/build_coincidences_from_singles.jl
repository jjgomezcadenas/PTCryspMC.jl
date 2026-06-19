#!/usr/bin/env julia
# Build the list-mode LOR file from a SINGLES stack (scripts/simulate_source_mt.jl output) —
# the production coincidence builder. One record per accepted coincidence = one LOR. The singles
# row already IS the formed per-gamma hit, so this fills `GammaAcc` directly (vs the full-stack
# accumulation in build_coincidences.jl) and runs the SAME selection core (src/coincidences.jl).
#
# Input: singles in either format (default prod/<tag>/singles.{h5,csv} by [output].format;
# --singles override; format by extension). Output: prod/<tag>/lors_{truth|det}.h5 — quantized
# Int16 columns, chunked + shuffle+deflate, with truth ∈ {true=0, scatter=1} (random=2 reserved
# for Step 5) and a `has_randoms=false` attribute (this is the LOR list MINUS randoms; the
# complete, time-ordered measurement = these ∪ the Step-5 random LORs, merged once times exist).
#
# Detector response comes from [detector] (all OFF → truth mode), identical to build_coincidences.
#
# Run from the repo root:
#   julia --project=. scripts/build_coincidences_from_singles.jl --config runs/sphere_water_csi.toml

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
    n_pair::Int; n_true::Int
end
LorState() = LorState(-1, 0, 0.0, 0.0, 0.0, 0, 0)

# Feed one singles row: close the previous event at a boundary (finish_event! → writer), then
# fill the gamma's accumulator directly. `g` = (GammaAcc, GammaAcc).
function feed_row!(st::LorState, g, w, resp, rng, ev::Int, gi::Int,
                   x::Float64, y::Float64, z::Float64, e::Float64, iz::Int, iphi::Int,
                   nblocks::Int, phscat::Bool, x0::Float64, y0::Float64, z0::Float64)
    if ev != st.cur_ev
        if st.cur_ev != -1
            emitted, is_true = finish_event!((a...) -> push_coincidence!(w, a...),
                                             st.cur_ev, g[1], g[2], st.ev_x0, st.ev_y0, st.ev_z0, resp, rng)
            st.n_pair += emitted; st.n_true += (emitted && is_true)
        end
        reset!(g[1]); reset!(g[2]); st.cur_ev = ev
    end
    ev > st.maxev && (st.maxev = ev)
    st.ev_x0 = x0; st.ev_y0 = y0; st.ev_z0 = z0
    (gi == 1 || gi == 2) || return
    fill_singles!(g[gi], x, y, z, e, iz, iphi, nblocks, phscat)
    return
end

function feed_singles_csv!(st, g, w, resp, rng, path)
    open(path, "r") do io
        header = split(strip(readline(io)), ',')
        col = Dict(String(h) => i for (i, h) in enumerate(header))
        for c in ("event_number", "gamma", "x_mm", "y_mm", "z_mm", "e_keV", "iz", "iphi",
                  "nblocks", "phantom_scatter", "x0_mm", "y0_mm", "z0_mm")
            haskey(col, c) || error("singles stack is missing column '$c'")
        end
        ie = col["event_number"]; ig = col["gamma"]; ix = col["x_mm"]; iy = col["y_mm"]; iz = col["z_mm"]
        iee = col["e_keV"]; iiz = col["iz"]; iip = col["iphi"]; inb = col["nblocks"]; iph = col["phantom_scatter"]
        ix0 = col["x0_mm"]; iy0 = col["y0_mm"]; iz0 = col["z0_mm"]
        for line in eachline(io)
            isempty(line) && continue
            f = split(line, ',')
            feed_row!(st, g, w, resp, rng, parse(Int, f[ie]), parse(Int, f[ig]),
                      parse(Float64, f[ix]), parse(Float64, f[iy]), parse(Float64, f[iz]),
                      parse(Float64, f[iee]), parse(Int, f[iiz]), parse(Int, f[iip]),
                      parse(Int, f[inb]), f[iph] == "1",
                      parse(Float64, f[ix0]), parse(Float64, f[iy0]), parse(Float64, f[iz0]))
        end
    end
end

function feed_singles_hdf5!(st, g, w, resp, rng, path)
    foreach_singles_hdf5(path) do b
        for i in 1:length(b)
            feed_row!(st, g, w, resp, rng, Int(b.event[i]), Int(b.gamma[i]),
                      decode_xyz(b.x[i]), decode_xyz(b.y[i]), decode_xyz(b.z[i]), decode_e(b.e[i]),
                      Int(b.iz[i]), Int(b.iphi[i]), Int(b.nblocks[i]), b.phantom_scatter[i] == 1,
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

    sigma_xyz = Float64(cfg_get(cfg, "detector", "sigma_xyz_mm", 0.0))
    eres      = Float64(cfg_get(cfg, "detector", "eres", 0.0))
    emin      = Float64(cfg_get(cfg, "detector", "emin_keV", 0.0))
    window    = Float64(cfg_get(cfg, "detector", "window_fwhm", 0.0))
    seed      = Int(cfg_get(cfg, "detector", "seed", 1234))
    window > 0.0 && eres <= 0.0 &&
        error("[detector].window_fwhm requires eres > 0 (the window width scales with the resolution)")
    resp = Response(sigma_xyz, eres, emin, window > 0.0, window * energy_fwhm(511.0, eres))
    mode = (sigma_xyz > 0.0 || eres > 0.0 || window > 0.0 || emin > 0.0) ? "det" : "truth"

    fmt = String(cfg_get(cfg, "output", "format", "csv"))
    default_in = joinpath(outdir, fmt == "hdf5" ? "singles.h5" : "singles.csv")
    singles = isempty(a["singles"]) ? default_in : rp(a["singles"])
    isfile(singles) || error("singles file '$singles' not found (run simulate_source_mt.jl first)")
    ishdf5 = endswith(singles, ".h5")
    out = isempty(a["out"]) ? joinpath(outdir, "lors_$mode.h5") : rp(a["out"])

    if mode == "det"
        ecut = window > 0.0 ? "window ±$(round(resp.win_half,digits=1)) keV" :
               (emin > 0.0 ? "Emin $(round(emin,digits=0)) keV" : "no energy cut")
        println("run '$tag' [det]: σ_xyz=$sigma_xyz mm, eres=$(round(100*eres,digits=1))% @511 keV, $ecut (seed $seed)")
    else
        println("run '$tag' [truth]: no smearing, no energy cut")
    end
    println("  singles ($(ishdf5 ? "hdf5" : "csv")): $singles")

    meta = Dict{String,Any}("scenario_tag" => tag, "mode" => mode, "has_randoms" => false,
        "sigma_xyz_mm" => sigma_xyz, "eres" => eres, "emin_keV" => emin, "window_fwhm" => window,
        "seed" => seed)
    w = CoincidenceWriter(out, meta)
    st = LorState()
    g = (GammaAcc(), GammaAcc())
    rng = MersenneTwister(seed)

    ishdf5 ? feed_singles_hdf5!(st, g, w, resp, rng, singles) :
             feed_singles_csv!(st, g, w, resp, rng, singles)
    if st.cur_ev != -1                                   # the final event
        emitted, is_true = finish_event!((a...) -> push_coincidence!(w, a...),
                                         st.cur_ev, g[1], g[2], st.ev_x0, st.ev_y0, st.ev_z0, resp, rng)
        st.n_pair += emitted; st.n_true += (emitted && is_true)
    end
    # Total annihilations for the acceptance denominator (the singles stack has no rows for
    # both-missed events). HDF5 carries the exact `nevents` attribute; CSV has no metadata, so
    # use the max event seen (≈ nevents). Stamp it onto the LOR file.
    nevents = ishdf5 ? Int(singles_hdf5_attr(singles, "nevents", 0)) : 0
    nevents > 0 || (nevents = st.maxev)
    set_lor_attr!(w, "nevents", nevents)
    nlor = close(w)

    println("wrote $nlor LORs -> $out")
    if nlor > 0
        denom = nevents > 0 ? nevents : st.n_pair
        tpct = round(100 * st.n_true / nlor, digits=1)
        println("  acceptance: $(round(100*nlor/max(denom,1), digits=1))% of $nevents annihilations")
        println("  truth split: $(st.n_true) true ($tpct%), $(nlor - st.n_true) scatter ($(round(100-tpct, digits=1))%)")
        println("  (truth ∈ {true,scatter}; random reserved for Step 5; has_randoms=false)")
    end
end

main()
