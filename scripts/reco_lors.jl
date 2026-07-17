#!/usr/bin/env julia
# Reco merge — the final detector-level list-mode file(s). Combines the truth coincidences
# (lors_truth.h5: true + scatter) and the randoms (randoms.h5), applying the SAME detector response
# + selection to both:
#   • smear each hit — position by σ_xyz, energy by FWHM(E) = a·√(511·E);
#   • energy-select the smeared energies (emin or a symmetric window);
#   • DT cut |t1−t2| ≤ τ (trims the trues' timing tail; randoms already satisfy it by construction);
#   • keep the truth flag (true=0 / scatter=1 / random=2) and the emitting isotope.
#
# GENERATION-2 (per-scenario timing): when [timing].t_ac_s + t_del_list are set, the transport spans
# the union window and this step BAND-CUTS each fixed acquisition scenario [t_del, t_del+t_ac] on the
# irradiation-end clock (t_decay_s), writing one prod/<tag>/lors_det_del<NNN>.h5 per scenario. Each
# hit is smeared ONCE and routed to every scenario whose window contains its decay (windows overlap;
# a decay is the same physical detection in each). Each shard is fully self-describing (dev/
# generation2_plan.md §5): generation="v2", the irradiation-end clock, the window, the embedded
# geometry, and the Mizuno washout g_i (NOT applied — carried for downstream). Absent → legacy single
# lors_det.h5 (t_decay_zero="acquisition_start").
#
#   julia --project=. scripts/reco_lors.jl --config runs/uniform_headep_ring1m_csi_v2.toml

using PTCryspMC
using ArgParse
using Random
import JSON

function parse_cli()
    s = ArgParseSettings(description="Reco merge: truth + randoms → detector-level lors_det.h5.")
    @add_arg_table! s begin
        "--config";  help = "run config TOML"; required = true
        "--truth";   help = "override the truth LOR path (default prod/<tag>/lors_truth.h5)"; default = ""
        "--randoms"; help = "override the randoms LOR path (default prod/<tag>/randoms.h5)"; default = ""
        "--out";     help = "override the output path (legacy single-file mode only)"; default = ""
    end
    parse_args(s)
end

# One output shard: its writer, its acquisition window [t1, t2] on the t_decay clock, its survivor
# tally (index = truth+1), and a label. Legacy mode uses a single sink with an infinite window.
struct RecoSink
    w::CoincidenceWriter
    t1::Float64
    t2::Float64
    counts::Vector{Int}
    label::String
end

# Stream one input LOR file through the detector pipeline once; route each survivor to every sink
# whose window contains its decay time. Returns rows read.
function reco_file!(sinks::Vector{RecoSink}, path, resp, tau::Float64, rng)::Int
    nread = 0
    foreach_coincidences_hdf5(path) do b
        for i in 1:length(b)
            nread += 1
            t1 = Float64(b.t1[i]); t2 = Float64(b.t2[i])
            tau > 0.0 && abs(t1 - t2) > tau && continue                      # DT cut
            if resp.sel_eff < 1.0        # model-2 selection tier: per-gamma Bernoulli (guarded
                k1 = rand(rng) < resp.sel_eff       # so sel_eff = 1 draws nothing — legacy
                k2 = rand(rng) < resp.sel_eff       # streams stay bit-identical)
                (k1 && k2) || continue
            end
            e1 = smear_energy(decode_e(b.e1[i]), resp.eres, rng)
            e2 = smear_energy(decode_e(b.e2[i]), resp.eres, rng)
            (pass_energy(e1, resp) && pass_energy(e2, resp)) || continue     # energy selection
            p1 = smear_hit((decode_xyz(b.x1[i]), decode_xyz(b.y1[i]), decode_xyz(b.z1[i])), resp, rng)
            p2 = smear_hit((decode_xyz(b.x2[i]), decode_xyz(b.y2[i]), decode_xyz(b.z2[i])), resp, rng)
            tr = Int(b.truth[i]); iso = Int(b.isotope[i]); td = Float64(b.t_decay[i])
            ev = Int(b.event[i]); tdc = b.t_decay[i]
            for s in sinks
                (td >= s.t1 && td <= s.t2) || continue                       # per-scenario band-cut
                push_coincidence!(s.w, ev, p1[1], p1[2], p1[3], e1, t1, Int(b.iz1[i]), Int(b.iphi1[i]), Int(b.nscat1[i]),
                    p2[1], p2[2], p2[3], e2, t2, Int(b.iz2[i]), Int(b.iphi2[i]), Int(b.nscat2[i]), Float64(b.dt[i]),
                    decode_xyz(b.x0[i]), decode_xyz(b.y0[i]), decode_xyz(b.z0[i]), tr, iso, tdc)
                s.counts[tr + 1] += 1
            end
        end
    end
    nread
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(prod_base(cfg)), tag)

    crystal = String(cfg_get(cfg, "transport", "crystal_material", "CsI"))
    mat = load_material(joinpath(REPO, "data"), crystal)

    sigma_cfg = Float64(cfg_get(cfg, "detector", "sigma_xyz_mm", NaN))
    eres_cfg  = Float64(cfg_get(cfg, "detector", "eres", NaN))
    sigma_xyz = isnan(sigma_cfg) ? mat.sigma_xyz_mm : sigma_cfg
    emin      = Float64(cfg_get(cfg, "detector", "reco_emin_keV", cfg_get(cfg, "detector", "emin_keV", 0.0)))
    window    = Float64(cfg_get(cfg, "detector", "window_fwhm", 0.0))
    seed      = Int(cfg_get(cfg, "detector", "seed", 1234))
    eres      = isnan(eres_cfg) ? mat.eres_a : eres_cfg
    window > 0.0 && eres <= 0.0 &&
        error("[detector].window_fwhm requires eres > 0 (window width scales with the resolution)")

    # Position model: 1 = single Gaussian σ_xyz (default), 2 = core/tail mixture (the crystal's
    # pos2 constants), with an optional selection tier ("none"/"80"/"60" — the per-gamma keep
    # efficiency 1.0/0.8/0.6, each with its own core fraction). Tiers are a model-2 concept.
    pos_model = Int(cfg_get(cfg, "detector", "pos_model", 1))
    seltier   = String(cfg_get(cfg, "detector", "selection", "none"))
    pos_model in (1, 2) || error("[detector].pos_model must be 1 or 2 (got $pos_model)")
    seltier in ("none", "80", "60") || error("[detector].selection must be \"none\", \"80\" or \"60\"")
    seltier != "none" && pos_model != 2 &&
        error("[detector].selection = \"$seltier\" requires pos_model = 2 (the tier sets the model-2 core fraction)")
    sel_eff, f_core = seltier == "none" ? (1.0, mat.pos2_f_none) :
                      seltier == "80"   ? (0.8, mat.pos2_f80)    : (0.6, mat.pos2_f60)
    pos_model == 2 && (mat.pos2_sigma1_mm <= 0.0 || f_core <= 0.0) &&
        error("pos_model = 2 needs pos2_sigma{1,2}_mm + pos2_f_core[\"$seltier\"] on crystal '$crystal' in materials.json")
    resp = Response(sigma_xyz, eres, emin, window > 0.0, window * energy_fwhm(511.0, eres),
                    pos_model, mat.pos2_sigma1_mm, mat.pos2_sigma2_mm, f_core, sel_eff)
    tau  = Float64(cfg_get(cfg, "timing", "tau_ns", 0.0))

    truth   = isempty(a["truth"])   ? joinpath(outdir, "lors_truth.h5") : rp(a["truth"])
    randoms = isempty(a["randoms"]) ? joinpath(outdir, "randoms.h5")    : rp(a["randoms"])
    isfile(truth)   || error("truth LORs '$truth' not found (run build_true_coincidences_from_singles.jl)")
    isfile(randoms) || error("randoms '$randoms' not found (run build_randoms_from_singles.jl)")
    tau > 0.0 || @warn "[timing].tau_ns ≤ 0: merging randoms with NO DT cut on the trues (inconsistent)"

    # Transport provenance (carried from the truth file so a shard alone regenerates the singles).
    nevents        = Int(singles_hdf5_attr(truth, "nevents", 0))
    transport_seed = Int(singles_hdf5_attr(truth, "transport_seed", -1))
    nchunks_tr     = Int(singles_hdf5_attr(truth, "nchunks", 0))

    # Scenario set (v2): fixed acquisitions [t_del, t_del+t_ac] on the irradiation-end clock. Absent
    # → legacy single lors_det.h5 (no band-cut).
    t_ac   = Float64(cfg_get(cfg, "timing", "t_ac_s", 0.0))
    t_dels = Float64[Float64(x) for x in cfg_get(cfg, "timing", "t_del_list", Float64[])]
    v2     = t_ac > 0.0 && !isempty(t_dels)

    # base attrs shared by every sink.
    base = Dict{String,Any}("scenario_tag" => tag, "mode" => "det", "has_randoms" => true,
        "crystal" => crystal, "seed" => seed, "t_relative_to_decay" => true,
        "sigma_xyz_mm" => sigma_xyz, "eres" => eres, "emin_keV" => emin, "window_fwhm" => window,
        "tau_ns" => tau, "pos_model" => pos_model)
    if pos_model == 2
        merge!(base, Dict{String,Any}("pos2_sigma1_mm" => mat.pos2_sigma1_mm,
            "pos2_sigma2_mm" => mat.pos2_sigma2_mm, "pos2_f_core" => f_core,
            "selection" => seltier, "sel_eff" => sel_eff))
    end

    ecut = window > 0.0 ? "window ±$(round(resp.win_half,digits=1)) keV" :
           (emin > 0.0 ? "Emin $(round(emin,digits=0)) keV" : "no energy cut")
    pos  = pos_model == 2 ?
        "pos2 σ1/σ2/f=$(mat.pos2_sigma1_mm)/$(mat.pos2_sigma2_mm)/$(round(f_core,digits=3)) mm" *
        (sel_eff < 1.0 ? " sel $(seltier)% (eff $sel_eff)" : "") :
        "σ_xyz=$sigma_xyz mm"
    println("run '$tag' reco: $pos, eres=$(round(100*eres,digits=1))%, $ecut, " *
            (tau > 0 ? "τ=$tau ns" : "no DT cut") * (v2 ? "  [v2: $(length(t_dels)) scenarios]" : ""))

    sinks = RecoSink[]
    if v2
        t_irr = Float64(singles_hdf5_attr(truth, "t_irr_s", 0.0))
        hls   = Vector{Float64}(singles_hdf5_attr(truth, "isotope_half_lives", Float64[]))
        isempty(hls) && error("v2 reco needs isotope_half_lives in $truth (API source)")
        # Embed the full geometry (self-containment) + key scalars.
        geomfile = rp(cfg_get(cfg, "geometry", "file", "geometry/geometry_head.json"))
        gjson    = read(geomfile, String)
        sc       = get(JSON.parse(gjson), "scanner", Dict{String,Any}())
        scname   = String(get(sc, "name", "scanner"))
        r_mm     = 10.0 * Float64(get(sc, "r_inner_cm", 0.0))
        afov_mm  = 20.0 * Float64(get(sc, "half_length_cm", 0.0))
        depth_mm = 10.0 * Float64(get(sc, "wall_thickness_cm", 0.0))
        for td in t_dels
            t1w = td; t2w = td + t_ac
            gvec = Float64[washout_gi(log(2.0) / h, t_irr, t1w, t2w) for h in hls]
            meta = merge(base, Dict{String,Any}(
                "generation" => "v2", "schema_version" => 2, "t_decay_zero" => "irradiation_end",
                "t_del_s" => td, "t_ac_s" => t_ac, "t1_s" => t1w, "t2_s" => t2w,
                "washout_model" => "mizuno_brain_3comp",
                "washout_fractions" => collect(MIZUNO_FRACTIONS), "washout_Thalf_s" => collect(MIZUNO_THALF_S),
                "washout_g" => gvec, "washout_applied" => false,
                "sigma_t_ns" => mat.sigma_t_ns, "eres_fwhm_511" => eres,
                "geometry_json" => gjson, "scanner_name" => scname,
                "ring_radius_mm" => r_mm, "afov_mm" => afov_mm, "crystal_depth_mm" => depth_mm))
            out = joinpath(outdir, "lors_det_del$(lpad(round(Int, td), 3, '0')).h5")
            push!(sinks, RecoSink(CoincidenceWriter(out, meta), t1w, t2w, zeros(Int, 3),
                                  "del$(round(Int, td))s [$(round(Int,t1w)),$(round(Int,t2w))]s"))
        end
    else
        out = isempty(a["out"]) ? joinpath(outdir, "lors_det.h5") : rp(a["out"])
        push!(sinks, RecoSink(CoincidenceWriter(out, base), -Inf, Inf, zeros(Int, 3), "lors_det"))
    end

    rng = MersenneTwister(seed)
    ntru = reco_file!(sinks, truth, resp, tau, rng)
    nrnd = reco_file!(sinks, randoms, resp, tau, rng)

    println("read $ntru truth + $nrnd random LORs")
    for s in sinks
        copy_provenance!(s.w, truth)                     # scenario/geometry provenance → self-describing
        set_lor_attr!(s.w, "nevents", nevents)
        transport_seed >= 0 && set_lor_attr!(s.w, "transport_seed", transport_seed)
        nchunks_tr > 0 && set_lor_attr!(s.w, "nchunks", nchunks_tr)
        nlor = close(s.w)
        c = s.counts
        pc(n) = nlor > 0 ? "$n ($(round(100*n/nlor, digits=1))%)" : "$n"
        println("  $(s.label): kept $nlor  →  $(pc(c[1])) true / $(pc(c[2])) scatter / $(pc(c[3])) random" *
                (nevents > 0 ? "  | acc $(round(100*nlor/nevents, digits=2))%" : ""))
    end
end

main()
