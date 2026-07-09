#!/usr/bin/env julia
# Reco merge — the final detector-level list-mode file. Combines the truth coincidences
# (lors_truth.h5: true + scatter) and the randoms (randoms.h5) into prod/<tag>/lors_det.h5,
# applying the SAME detector response + selection to both:
#   • smear each hit — position by σ_xyz, energy by FWHM(E) = a·√(511·E);
#   • energy-select the smeared energies (emin or a symmetric window);
#   • DT cut |t1−t2| ≤ τ (trims the trues' timing tail; randoms already satisfy it by construction);
#   • keep the truth flag (true=0 / scatter=1 / random=2).
# The output is the union, flagged, NOT time-sorted: the LOR list is order-independent for
# reconstruction, and any consumer that wants a chronological stream sorts by the recomputable
# absolute time event_time(event)·1e9 + min(t1,t2). Streaming (O(1) memory).
#
# eres defaults to the crystal's own eres_a unless [detector].eres is set; τ from [timing].tau_ns.
#
# Note: trues and randoms are not mutually exclusive (a single may appear in both a true and a
# random LOR) — the multiples approximation, negligible at our rates; see src/randoms.jl.
#
#   julia --project=. scripts/reco_lors.jl --config runs/sphere_water_csi.toml

using PTCryspMC
using ArgParse
using Random

function parse_cli()
    s = ArgParseSettings(description="Reco merge: truth + randoms → detector-level lors_det.h5.")
    @add_arg_table! s begin
        "--config";  help = "run config TOML"; required = true
        "--truth";   help = "override the truth LOR path (default prod/<tag>/lors_truth.h5)"; default = ""
        "--randoms"; help = "override the randoms LOR path (default prod/<tag>/randoms.h5)"; default = ""
        "--out";     help = "override the output path (default prod/<tag>/lors_det.h5)"; default = ""
    end
    parse_args(s)
end

# Stream one input LOR file through the detector pipeline, appending survivors to the writer.
# `counts` tallies survivors by truth code (index = truth+1). Returns rows read.
function reco_file!(w, path, resp, tau::Float64, rng, counts::Vector{Int})::Int
    nread = 0
    foreach_coincidences_hdf5(path) do b
        for i in 1:length(b)
            nread += 1
            t1 = Float64(b.t1[i]); t2 = Float64(b.t2[i])
            tau > 0.0 && abs(t1 - t2) > tau && continue                      # DT cut
            e1 = smear_energy(decode_e(b.e1[i]), resp.eres, rng)
            e2 = smear_energy(decode_e(b.e2[i]), resp.eres, rng)
            (pass_energy(e1, resp) && pass_energy(e2, resp)) || continue     # energy selection
            x1, y1, z1 = smear_position((decode_xyz(b.x1[i]), decode_xyz(b.y1[i]), decode_xyz(b.z1[i])), resp.sigma_xyz, rng)
            x2, y2, z2 = smear_position((decode_xyz(b.x2[i]), decode_xyz(b.y2[i]), decode_xyz(b.z2[i])), resp.sigma_xyz, rng)
            tr = Int(b.truth[i])
            push_coincidence!(w, Int(b.event[i]), x1, y1, z1, e1, t1, Int(b.iz1[i]), Int(b.iphi1[i]), Int(b.nscat1[i]),
                x2, y2, z2, e2, t2, Int(b.iz2[i]), Int(b.iphi2[i]), Int(b.nscat2[i]), Float64(b.dt[i]),
                decode_xyz(b.x0[i]), decode_xyz(b.y0[i]), decode_xyz(b.z0[i]), tr, b.t_decay[i])
            counts[tr + 1] += 1
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

    sigma_cfg = Float64(cfg_get(cfg, "detector", "sigma_xyz_mm", NaN)) # NaN = absent → crystal default
    eres_cfg  = Float64(cfg_get(cfg, "detector", "eres", NaN))         # NaN = absent → crystal default
    sigma_xyz = isnan(sigma_cfg) ? mat.sigma_xyz_mm : sigma_cfg        # per-crystal σ (CRYSP 3.5 mm FWHM)
    # The reconstructed LORs use a HIGHER lower energy cut than the spectrum studies (which keep
    # emin_keV low to see the Compton shoulder): reco_emin_keV selects the photopeak region (~a few σ
    # below 511 for the worst-resolution crystal) and shrinks the file. Falls back to emin_keV if unset.
    emin      = Float64(cfg_get(cfg, "detector", "reco_emin_keV", cfg_get(cfg, "detector", "emin_keV", 0.0)))
    window    = Float64(cfg_get(cfg, "detector", "window_fwhm", 0.0))
    seed      = Int(cfg_get(cfg, "detector", "seed", 1234))
    eres      = isnan(eres_cfg) ? mat.eres_a : eres_cfg                # name the crystal → its eres
    window > 0.0 && eres <= 0.0 &&
        error("[detector].window_fwhm requires eres > 0 (window width scales with the resolution)")
    resp = Response(sigma_xyz, eres, emin, window > 0.0, window * energy_fwhm(511.0, eres))
    tau  = Float64(cfg_get(cfg, "timing", "tau_ns", 0.0))

    truth   = isempty(a["truth"])   ? joinpath(outdir, "lors_truth.h5") : rp(a["truth"])
    randoms = isempty(a["randoms"]) ? joinpath(outdir, "randoms.h5")    : rp(a["randoms"])
    out     = isempty(a["out"])     ? joinpath(outdir, "lors_det.h5")   : rp(a["out"])
    isfile(truth)   || error("truth LORs '$truth' not found (run build_true_coincidences_from_singles.jl)")
    isfile(randoms) || error("randoms '$randoms' not found (run build_randoms_from_singles.jl)")
    # Merging randoms (paired within some τ) but applying no DT cut to the trues would be an
    # inconsistent lors_det — guard against a config that cleared [timing].tau_ns after the randoms
    # were built. (build_randoms itself errors on τ ≤ 0.)
    tau > 0.0 || @warn "[timing].tau_ns ≤ 0: merging randoms with NO DT cut on the trues (inconsistent)"

    ecut = window > 0.0 ? "window ±$(round(resp.win_half,digits=1)) keV" :
           (emin > 0.0 ? "Emin $(round(emin,digits=0)) keV" : "no energy cut")
    eresrc = isnan(eres_cfg) ? " (from $crystal)" : ""
    println("run '$tag' reco: σ_xyz=$sigma_xyz mm, eres=$(round(100*eres,digits=1))%$eresrc, $ecut, " *
            (tau > 0 ? "τ=$tau ns" : "no DT cut") * " (seed $seed)")

    nevents = Int(singles_hdf5_attr(truth, "nevents", 0))
    # The transport recipe (seed + chunk count), carried from lors_truth so lors_det.h5 alone
    # suffices to regenerate the pruned singles. Absent (-1/0) when the truth came from CSV singles.
    transport_seed = Int(singles_hdf5_attr(truth, "transport_seed", -1))
    nchunks_tr     = Int(singles_hdf5_attr(truth, "nchunks", 0))
    meta = Dict{String,Any}("scenario_tag" => tag, "mode" => "det", "has_randoms" => true,
        "crystal" => crystal, "seed" => seed, "t_relative_to_decay" => true,
        "sigma_xyz_mm" => sigma_xyz, "eres" => eres, "emin_keV" => emin, "window_fwhm" => window,
        "tau_ns" => tau)
    w = CoincidenceWriter(out, meta)
    rng = MersenneTwister(seed)
    counts = zeros(Int, 3)                                # [true, scatter, random]

    ntru = reco_file!(w, truth, resp, tau, rng, counts)
    nrnd = reco_file!(w, randoms, resp, tau, rng, counts)
    set_lor_attr!(w, "nevents", nevents)
    transport_seed >= 0 && set_lor_attr!(w, "transport_seed", transport_seed)
    nchunks_tr > 0 && set_lor_attr!(w, "nchunks", nchunks_tr)
    copy_provenance!(w, truth)           # carry the source/geometry provenance → lors_det is self-describing
    nlor = close(w)

    println("read $ntru truth + $nrnd random LORs → kept $nlor → $out")
    if nlor > 0
        pc(n) = "$n ($(round(100*n/nlor, digits=1))%)"
        println("  truth split:  $(pc(counts[1])) true / $(pc(counts[2])) scatter / $(pc(counts[3])) random")
        println("  randoms / (true+scatter): $(round(100*counts[3]/max(counts[1]+counts[2],1), digits=2))%")
        nevents > 0 && println("  acceptance:   $(round(100*nlor/nevents, digits=2))% of $nevents annihilations")
    end
end

main()
