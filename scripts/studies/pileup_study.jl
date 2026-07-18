#!/usr/bin/env julia
# Pileup (front-end integration) study — quantifies the "one annihilation at a time" assumption.
#
# The front-end integrates each block's signal for T_int after a deposit (CsI 3 µs; BGO 4 µs at
# 195 K, 10 µs at 77 K — a few times the scintillation decay). A hit is PILED UP if another gamma
# deposits in the SAME block while its gate is open; a LOR is piled up if either of its two hits
# is. With a per-block singles rate r (Poisson arrivals) the per-LOR probability is
#
#     P_LOR(t) = 1 − exp(−2 · r(t) · T_int),    r(t) = A(t) · Y · s_block
#
# where A(t) is the scenario activity at time t after irradiation end (analytic, from the decay
# budget: A(t) = Σ_j λ_j N⁰_j e^{−λ_j t}), Y the detected singles per decay, and s_block the
# block's share of the singles (mean = 1/nblocks; hottest block measured). Y and s_block are
# MEASURED here by running the actual transport (the same `singles_chunk!` the production driver
# uses) at a small dose (DOSE_MEAS, ~1e6 decays) and counting every deposit above the 10 keV
# transport cutoff per block — every deposit corrupts an open gate, so the full singles rate
# (not the photopeak or trigger rate) is the conservative, correct rate.
#
# Outputs (studies/pileup/):
#   pileup_table.csv    one row per scanner × crystal readout (Y, s_hot, rates, P_LOR at the
#                       acquisition starts and LOR-weighted over each 300 s window)
#   activity_curve.csv  A(t) per isotope + total, t = 0…600 s (for py/fig_pileup.py)
# and prints the human table + the LaTeX table body for latex/dead_time.tex.
#
#   julia --project=. scripts/studies/pileup_study.jl

using PTCryspMC
using Random
using Printf

const REPO      = normpath(joinpath(@__DIR__, "..", ".."))
const SCND      = joinpath(REPO, "..", "ptcrysp-scenarios", "scenarios", "uniform_headep_sobp_1e8")
const OUTDIR    = joinpath(REPO, "studies", "pileup")
const BUDGET    = "fast"
const DOSE      = 1.0      # Gy — the activity the rates are quoted at (dose-linear)
const DOSE_MEAS = 0.02     # Gy — the measurement runs (~1e6 decays: ±few % on the hot block)
const SEED      = 1234
const NCHUNKS   = 8
const E0        = 0.511    # MeV
const CUT_MEV   = 0.010    # transport cutoff — every deposit above it counts
const ACOL      = 0.5      # deg FWHM
const T_DEL     = [120.0, 180.0, 300.0]   # the v2 acquisition starts (irradiation-end clock)
const T_AC      = 300.0

# Front-end integration time per crystal readout [µs] (electronics: a few × scintillation decay).
const T_INT_US = Dict("CsI" => 3.0, "BGO_195K" => 4.0, "BGO_77K" => 10.0)

# Scanners: geometry file, short label, crystal readouts. The two BGO temperatures share the
# attenuation (same xcom_BGO table), so one transport serves both rows — only T_int differs.
const SCANNERS = [
    ("geometry_head_csi_2x0.json",      "ring 1 m",     ["CsI"]),
    ("geometry_ring1m_bgo_2x0.json",    "ring 1 m",     ["BGO_195K", "BGO_77K"]),
    ("geometry_r35_50cm_csi_2x0.json",  "r35 50 cm",    ["CsI"]),
    ("geometry_r40_50cm_bgo_2x0.json",  "r40 50 cm",    ["BGO_195K", "BGO_77K"]),
    ("geometry_r35_35cm_csi_2x0.json",  "r35 35 cm",    ["CsI"]),
    ("geometry_r40_35cm_bgo_2x0.json",  "r40 35 cm",    ["BGO_195K", "BGO_77K"]),
    ("geometry_chs_csi_2x0.json",       "CHS",          ["CsI"]),
    ("geometry_chs_bgo_2x0.json",       "CHS",          ["BGO_195K", "BGO_77K"]),
]

# ---------------------------------------------------------------- activity (analytic)

"Tiny CSV reader: header names → column vectors of String."
function read_csv_cols(path::AbstractString)
    lines = readlines(path)
    hdr   = split(strip(lines[1]), ",")
    cols  = Dict(String(h) => String[] for h in hdr)
    for ln in lines[2:end]
        isempty(strip(ln)) && continue
        for (h, v) in zip(hdr, split(strip(ln), ","))
            push!(cols[String(h)], String(v))
        end
    end
    cols
end

"""
Back the irradiation-end populations N⁰_j out of the sampling budget (decays in
[t_del_b, t_del_b + t_meas], at the budget's reference dose) and return
(names, λ [1/s], N⁰ at `dose` Gy). A(t) = Σ λ_j N⁰_j exp(−λ_j t).
"""
function activity_terms(scndir::AbstractString, budget::AbstractString, dose::Float64)
    iso   = read_csv_cols(joinpath(scndir, "isotopes.csv"))
    bud   = read_csv_cols(joinpath(scndir, "sampling_budget_$(budget).csv"))
    meta  = read_csv_cols(joinpath(scndir, "sampling_budget_$(budget)_meta.csv"))
    t_meas  = parse(Float64, meta["t_meas_s"][1])
    t_del_b = haskey(meta, "t_del_s") ? parse(Float64, meta["t_del_s"][1]) : 0.0
    refdose = parse(Float64, meta["dose_Gy"][1])
    thalf = Dict(id => parse(Float64, th) for (id, th) in zip(iso["isotope_id"], iso["half_life_s"]))
    names = Dict(id => nm for (id, nm) in zip(iso["isotope_id"], iso["name"]))
    nm = String[]; lam = Float64[]; n0 = Float64[]
    for (id, ne) in zip(bud["isotope_id"], bud["N_expected"])
        λ = log(2.0) / thalf[id]
        f = exp(-λ * t_del_b) - exp(-λ * (t_del_b + t_meas))   # budget-window factor
        push!(nm, names[id]); push!(lam, λ)
        push!(n0, parse(Float64, ne) * (dose / refdose) / f)
    end
    nm, lam, n0
end

activity(lam, n0, t) = sum(λ * N * exp(-λ * t) for (λ, N) in zip(lam, n0))   # Bq

"LOR-weighted effective activity over [t0, t0+tac]: ∫A²dt / ∫A dt (pileup ∝ A, LOR rate ∝ A)."
function activity_eff(lam, n0, t0, tac; dt=0.5)
    ts = t0:dt:(t0 + tac)
    sum(activity(lam, n0, t)^2 for t in ts) / sum(activity(lam, n0, t) for t in ts)
end

# ---------------------------------------------------------------- rate measurement

"Transport ~DOSE_MEAS worth of decays through `geomfile` with `crystal`; per-block hit counts."
function measure_blocks(geomfile::String, crystal::String, mats, scn, src)
    geom0 = load_geometry(joinpath(REPO, "geometry", geomfile), mats)
    sc0   = geom0.scanner
    sc    = Scanner(PhysicalVolume(LogicalVolume(name(sc0.volume), solid(sc0.volume),
                                                 mats[crystal]), sc0.volume.position),
                    sc0.n_phi, sc0.n_z)
    geom  = Geometry(geom0.world, scn.phantom, sc)
    N     = length(src.points)
    counts = Dict{Tuple{Int,Int},Int}()
    for (c, rg) in enumerate(chunk_ranges(N, NCHUNKS))
        rng = MersenneTwister(SEED + c - 1)
        singles_chunk!(geom, src, E0, CUT_MEV, ACOL, rg, rng, mats[crystal]) do ev, g, s, pos0, t_rel, iso
            counts[(s.iz, s.iphi)] = get(counts, (s.iz, s.iphi), 0) + 1
        end
    end
    (nblocks = nblocks(sc), n_decays = N, total = sum(values(counts)),
     hot = maximum(values(counts)), scanner = name(sc.volume))
end

# ---------------------------------------------------------------- run

p_lor(r, T_us) = 1.0 - exp(-2.0 * r * T_us * 1e-6)

function main()
    mkpath(OUTDIR)
    mats = load_materials(joinpath(REPO, "data"))
    nm, lam, n0 = activity_terms(SCND, BUDGET, DOSE)
    a(t) = activity(lam, n0, t)
    @printf("scenario activity at %.1f Gy: %.1f kBq at irradiation end", DOSE, a(0.0) / 1e3)
    for td in T_DEL
        @printf(", %.1f kBq at t=%.0f s", a(td) / 1e3, td)
    end
    println(".")

    # A(t) curve for the figure.
    open(joinpath(OUTDIR, "activity_curve.csv"), "w") do io
        println(io, "t_s," * join(["A_$(x)_kBq" for x in nm], ",") * ",A_total_kBq")
        for t in 0.0:2.0:600.0
            terms = [λ * N * exp(-λ * t) / 1e3 for (λ, N) in zip(lam, n0)]
            println(io, join(vcat([t], round.(terms, digits=4), [round(sum(terms), digits=4)]), ","))
        end
    end

    # One measurement source, shared by every scanner (the survey convention).
    scn = load_scenario(SCND, mats; budget=BUDGET, dose_Gy=DOSE_MEAS, keep_escaped=false,
                        center_on="tumour", t_window=(minimum(T_DEL), maximum(T_DEL) + T_AC))
    src = materialize_api_source(scn; master_seed=1, realization=0)

    rows = NamedTuple[]
    for (gf, label, crystals) in SCANNERS
        m = measure_blocks(gf, crystals[1], mats, scn, src)
        Y, s_hot = m.total / m.n_decays, m.hot / m.total
        @printf("%-34s %7d decays -> %8d singles (Y=%.3f), %d blocks, hot share %.3f%%\n",
                m.scanner, m.n_decays, m.total, Y, m.nblocks, 100 * s_hot)
        for cr in crystals
            T = T_INT_US[cr]
            r_mean(t) = a(t) * Y / m.nblocks
            r_hot(t)  = a(t) * Y * s_hot
            push!(rows, (scanner = m.scanner, label = label, crystal = cr, T_int_us = T,
                nblocks = m.nblocks, singles_per_decay = round(Y, digits=4),
                hot_frac = round(s_hot, digits=6),
                S_kHz_t0 = round(a(0.0) * Y / 1e3, digits=1),
                P_hot_t0 = p_lor(r_hot(0.0), T),
                (Symbol("P_mean_del$(Int(td))") => p_lor(r_mean(td), T) for td in T_DEL)...,
                (Symbol("P_hot_del$(Int(td))") => p_lor(r_hot(td), T) for td in T_DEL)...,
                (Symbol("P_hot_avg$(Int(td))") =>
                    p_lor(activity_eff(lam, n0, td, T_AC) * Y * s_hot, T) for td in T_DEL)...))
        end
    end

    open(joinpath(OUTDIR, "pileup_table.csv"), "w") do io
        ks = keys(rows[1])
        println(io, join(ks, ","))
        for r in rows
            println(io, join([r[k] isa Float64 && startswith(String(k), "P_") ?
                              @sprintf("%.6f", r[k]) : r[k] for k in ks], ","))
        end
    end

    println("\n--- P_LOR [%]: mean/hot block at the t_del=120 s acquisition start; " *
            "hot at irradiation end; hot LOR-weighted over [120,420] s")
    @printf("%-12s %-9s %5s %7s | %10s %9s %9s %9s\n",
            "scanner", "crystal", "T_us", "blocks", "mean@120", "hot@120", "hot@0", "hot avg120")
    for r in rows
        @printf("%-12s %-9s %5.1f %7d | %9.2f%% %8.2f%% %8.2f%% %8.2f%%\n",
                r.label, r.crystal, r.T_int_us, r.nblocks, 100 * r.P_mean_del120,
                100 * r.P_hot_del120, 100 * r.P_hot_t0, 100 * r.P_hot_avg120)
    end

    println("\n--- LaTeX table body (latex/dead_time.tex)")
    for r in rows
        @printf("%s & %s & %.0f & %d & %.2f & %.2f & %.2f\\%% & %.2f\\%% & %.2f\\%% & %.2f\\%% \\\\\n",
                r.label, replace(r.crystal, "_" => "\\_"), r.T_int_us, r.nblocks,
                r.singles_per_decay, 100 * r.hot_frac, 100 * r.P_mean_del120,
                100 * r.P_hot_del120, 100 * r.P_hot_t0, 100 * r.P_hot_avg120)
    end
end

main()
