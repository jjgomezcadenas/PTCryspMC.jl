#!/usr/bin/env julia
# CTR study — measure the coincidence timing of a scanner from MC and derive the window τ.
#
# Streams a truth LOR file (build_true_coincidences_from_singles.jl output; NO DT cut applied
# there, so the distribution is unbiased) and characterizes, per LOR, the two signed timing
# quantities:
#   • Δt_raw = t1 − t2              — what a DAQ windows on: detector jitter ⊕ the geometric TOF
#                                     spread of the source across the phantom. τ is set on THIS.
#   • dt     = (t1 − t2) − TOF_diff — the TOF-corrected residual: the detector-only CTR
#                                     (crystal characterization; undefined for randoms).
# For each it reports mean, σ (std), a robust σ (half-width of the central 68.27%), and the
# FWHM (histogram-based; CTR by the usual convention), plus the true-coincidence containment
# vs τ at σ-multiples and on an absolute grid. The first-photon jitter is exponential, so Δt is
# Laplace-like — heavier-tailed than a Gaussian — and the containment column, not the σ multiple,
# is what the τ choice should read.
#
# Output: a summary table (stdout), studies/ctr/<tag>_ctr.csv (the containment table) and
# studies/ctr/<tag>_ctr_hist.csv (binned Δt_raw and dt histograms, for plotting).
#
#   julia --project=. scripts/studies/ctr_study.jl --config runs/uniform_headep_bgo77k_api.toml

using PTCryspMC
using ArgParse
using Printf

function parse_cli()
    s = ArgParseSettings(description="Measure CTR + containment vs τ from a truth LOR file.")
    @add_arg_table! s begin
        "--config"; help = "run config TOML"; required = true
        "--lors";   help = "override the LOR path (default prod/<tag>/lors_truth.h5)"; default = ""
    end
    parse_args(s)
end

"Robust σ: half-width of the central 68.27% of the SORTED sample (equals σ for a Gaussian)."
robust_sigma(sorted::Vector{Float32}) = begin
    n = length(sorted)
    lo = sorted[clamp(round(Int, 0.15865 * n), 1, n)]
    hi = sorted[clamp(round(Int, 0.84135 * n), 1, n)]
    (hi - lo) / 2
end

"FWHM from a fine histogram of `v` (already sorted): width of the region above half the peak count."
function fwhm_hist(sorted::Vector{Float32}; nbins::Int=4000)
    n = length(sorted)
    lo = sorted[clamp(round(Int, 0.001 * n), 1, n)]     # clip the extreme tails for the binning range
    hi = sorted[clamp(round(Int, 0.999 * n), 1, n)]
    hi > lo || return 0.0
    Δ = (hi - lo) / nbins
    counts = zeros(Int, nbins)
    for x in sorted
        lo <= x <= hi && (counts[clamp(floor(Int, (x - lo) / Δ) + 1, 1, nbins)] += 1)
    end
    cmax = maximum(counts); half = cmax / 2
    i1 = findfirst(>=(half), counts); i2 = findlast(>=(half), counts)
    (i2 - i1 + 1) * Δ
end

"Containment: fraction of the SORTED |v| sample ≤ τ (binary search)."
containment(sorted_abs::Vector{Float32}, τ) = searchsortedlast(sorted_abs, Float32(τ)) / length(sorted_abs)

function summarize(name, v::Vector{Float32})
    sort!(v)
    n  = length(v)
    μ  = sum(Float64.(v)) / n
    σ  = sqrt(sum(x -> (Float64(x) - μ)^2, v) / (n - 1))
    σr = robust_sigma(v)
    fw = fwhm_hist(v)
    @printf("  %-28s mean %+.4f   σ %.4f   σ_robust %.4f   FWHM %.4f  [ns]\n", name, μ, σ, σr, fw)
    (μ=μ, σ=σ, σr=σr, fwhm=fw)
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    tag = run_tag(cfg, a["config"])
    lors = isempty(a["lors"]) ? joinpath(rp(prod_base(cfg)), tag, "lors_truth.h5") : rp(a["lors"])
    isfile(lors) || error("LOR file '$lors' not found (run the chain to lors_truth first)")
    outdir = joinpath(REPO, "studies", "ctr"); mkpath(outdir)
    crystal = String(cfg_get(cfg, "transport", "crystal_material", "?"))

    # Condition on the reco energy selection: the jitter scales as 1/E_dep (N_det = yield·E·pde),
    # so the uncut truth file contains few-keV hits with jitter of 10²–10⁵ ns that the energy cut
    # removes before any DT window applies. The window population is the energy-selected one.
    emin = Float64(cfg_get(cfg, "detector", "reco_emin_keV", cfg_get(cfg, "detector", "emin_keV", 0.0)))

    draw = Float32[]; dres = Float32[]
    ntot = 0
    foreach_coincidences_hdf5(lors) do b
        for i in 1:length(b)
            ntot += 1
            (decode_e(b.e1[i]) >= emin && decode_e(b.e2[i]) >= emin) || continue
            push!(draw, b.t1[i] - b.t2[i])   # signed raw Δt (the window variable)
            push!(dres, b.dt[i])             # signed TOF-corrected residual (the CTR variable)
        end
    end
    n = length(draw)
    println("CTR study '$tag' (crystal $crystal): $n of $ntot truth LORs pass E ≥ $emin keV (both gammas)")

    r  = summarize("Δt_raw = t1−t2 (window var)", draw)
    c  = summarize("dt = Δt−TOF_diff (CTR)", dres)

    # Containment vs τ: multiples of the ROBUST raw σ (the std is tail-inflated — the jitter is
    # exponential, so Δt is Laplace-like), then an absolute grid.
    abs_draw = sort!(abs.(draw))
    println("  containment of trues, |t1−t2| ≤ τ  (σ = robust raw σ = $(round(r.σr, digits=4)) ns):")
    rows = Tuple{String,Float64,Float64}[]
    for k in (1.0, 2.0, 2.5, 3.0, 4.0, 5.0)
        τ = k * r.σr
        push!(rows, (@sprintf("%.1fσ", k), τ, containment(abs_draw, τ)))
    end
    for τ in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 10.0)
        push!(rows, ("", Float64(τ), containment(abs_draw, τ)))
    end
    for (lbl, τ, frac) in rows
        @printf("    %-5s τ = %6.3f ns   keeps %8.4f %%   (loses %.4f %%)\n", lbl, τ, 100frac, 100(1 - frac))
    end

    # CSV: the containment table + the summary moments as header comments.
    csv = joinpath(outdir, "$(tag)_ctr.csv")
    open(csv, "w") do io
        @printf(io, "# crystal=%s n=%d\n", crystal, n)
        @printf(io, "# raw: mean=%.6f sigma=%.6f sigma_robust=%.6f fwhm=%.6f\n", r.μ, r.σ, r.σr, r.fwhm)
        @printf(io, "# ctr: mean=%.6f sigma=%.6f sigma_robust=%.6f fwhm=%.6f\n", c.μ, c.σ, c.σr, c.fwhm)
        println(io, "label,tau_ns,containment")
        for (lbl, τ, frac) in rows
            @printf(io, "%s,%.6f,%.8f\n", lbl, τ, frac)
        end
    end

    # Histograms for plotting: 10 ps bins over ±6σ of each variable.
    hist = joinpath(outdir, "$(tag)_ctr_hist.csv")
    open(hist, "w") do io
        println(io, "var,bin_center_ns,count")
        for (name, v, σ) in (("raw", draw, r.σ), ("ctr", dres, c.σ))
            lim = 6σ; nb = max(200, min(4000, round(Int, 2lim / 0.01)))
            Δ = 2lim / nb; counts = zeros(Int, nb)
            for x in v
                -lim <= x < lim && (counts[floor(Int, (x + lim) / Δ) + 1] += 1)
            end
            for i in 1:nb
                counts[i] > 0 && @printf(io, "%s,%.5f,%d\n", name, -lim + (i - 0.5) * Δ, counts[i])
            end
        end
    end
    println("  wrote $csv")
    println("  wrote $hist")
    println("  → choose τ = 2σ or 3σ from the containment column (Laplace tails: σ-multiples read low).")
end

main()
