#!/usr/bin/env julia
# QA for the API scenario reader: load a frozen ptcryspg4 scenario with `load_scenario` and report
# what the detector sim will see — the per-isotope annihilation-point pools, the expected decay
# budget N_expected (rescaled to the requested dose), the kept fraction f_inside and the escaped
# positrons dropped, plus the effective detectable count M_j = N_expected · f_inside the source
# (step 5) will Poisson-draw, and the per-isotope decay timing in the measurement window (the
# truncated exponential the randoms pass samples). Runs sanity gates (fractions in range, every
# budget isotope has a pool, the escaped fraction is small, each isotope's sampled decay-time mean
# matches the analytic truncated-exponential) and errors if any is violated, so it is a genuine check.
#
#   julia --project=. scripts/tests/check_scenario.jl --scenario ~/Projects/ptcrysp-scenarios/scenarios/uniform_headep_sobp_1e8
#   julia --project=. scripts/tests/check_scenario.jl --scenario <dir> --budget offline --dose 2.0

using PTCryspMC
using ArgParse
using Printf

function parse_cli()
    s = ArgParseSettings(description="Check the API scenario reader against a frozen scenario.")
    @add_arg_table! s begin
        "--scenario";     help = "scenario directory (holds emitters.csv, phantom_regions.csv, …)"; required = true
        "--budget";       help = "timing budget: fast | inroom | offline"; default = "fast"
        "--dose";         help = "clinical dose [Gy] (linear rescale of N_expected)"; arg_type = Float64; default = 1.0
        "--keep-escaped"; help = "keep escaped positrons (outside the phantom) instead of dropping"; action = :store_true
        "--data";         help = "materials/XCOM directory"; default = "data"
    end
    parse_args(s)
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    isdir(rp(a["scenario"])) || error("scenario directory '$(a["scenario"])' not found")
    mats = load_materials(rp(a["data"]))

    t = @elapsed s = load_scenario(rp(a["scenario"]), mats;
                                   budget=a["budget"], dose_Gy=a["dose"], keep_escaped=a["keep-escaped"])

    @printf("scenario '%s'  budget=%s  dose=%.3g Gy  (loaded in %.1fs)\n",
            s.name, s.budget, s.dose_Gy, t)
    @printf("  phantom: %s %s, material %s | t_meas=%.0f s | keep_escaped=%s\n",
            name(s.phantom), nameof(typeof(solid(s.phantom))), material(s.phantom).name,
            s.t_meas_s, s.provenance["keep_escaped"])
    @printf("  provenance: geometry=%s  g4=%s  physics=%s  n_protons=%s\n",
            s.provenance["geometry"], s.provenance["geant4_version"],
            s.provenance["physics_list"], s.provenance["n_protons"])

    n = length(s.pools)
    @printf("\n  %-5s %12s %12s %8s %9s %12s\n",
            "iso", "pool", "N_exp", "f_in", "dropped", "M_j=N·f_in")
    tot_exp = 0.0; tot_eff = 0.0; tot_pool = 0; tot_drop = 0
    for j in 1:n
        eff = s.n_expected[j] * s.f_inside[j]
        @printf("  %-5s %12d %12.4e %8.4f %9d %12.4e\n",
                s.isotopes[j].name, length(s.pools[j]), s.n_expected[j],
                s.f_inside[j], s.n_dropped[j], eff)
        tot_exp += s.n_expected[j]; tot_eff += eff
        tot_pool += length(s.pools[j]); tot_drop += s.n_dropped[j]
    end
    n_emit = tot_pool + tot_drop
    @printf("  %-5s %12d %12.4e %8s %9d %12.4e\n",
            "all", tot_pool, tot_exp, "", tot_drop, tot_eff)
    @printf("\n  emitters read: %d  |  escaped dropped: %d (%.3f%%)\n",
            n_emit, tot_drop, n_emit > 0 ? 100 * tot_drop / n_emit : 0.0)
    @printf("  total N_expected = %.4e  |  effective detectable ΣM_j = %.4e\n", tot_exp, tot_eff)

    # Per-isotope decay timing in the measurement window (the truncated-exponential the randoms
    # pass samples): short-lived isotopes front-load (mean → 1/λ), long-lived flatten (mean → T/2).
    models = scenario_activity_models(s; seed=1234)
    T = s.t_meas_s; ndraw = 50_000
    @printf("\n  per-isotope decay time in the %.0f s window (%d draws):\n", T, ndraw)
    @printf("  %-5s %10s %10s %12s %12s %8s\n", "iso", "T_half[s]", "1/λ[s]", "mean_samp", "mean_exact", "samp/ex")
    bad_t = String[]
    for j in 1:length(models)
        λ = models[j].λ
        acc = 0.0
        for ev in 1:ndraw; acc += event_time(models[j], ev); end
        samp  = acc / ndraw
        exact = 1.0 / λ - T * exp(-λ * T) / (1.0 - exp(-λ * T))   # truncated-exp mean on [0,T]
        @printf("  %-5s %10.1f %10.1f %12.1f %12.1f %8.4f\n",
                s.isotopes[j].name, s.isotopes[j].half_life_s, 1.0 / λ, samp, exact, samp / exact)
        (isapprox(samp, exact; rtol=0.05) && 0.0 <= samp <= T) || push!(bad_t, s.isotopes[j].name)
    end

    # ---- sanity gates (error on violation) ----------------------------------
    isempty(bad_t) ||
        error("isotopes whose sampled decay-time mean misses the truncated-exponential: $(join(bad_t, ", "))")
    for j in 1:n
        0.0 < s.f_inside[j] <= 1.0 ||
            error("isotope $(s.isotopes[j].name): f_inside $(s.f_inside[j]) out of (0,1]")
        s.n_expected[j] >= 0.0 || error("isotope $(s.isotopes[j].name): negative N_expected")
        (length(s.pools[j]) > 0 || s.n_expected[j] == 0.0) ||
            error("isotope $(s.isotopes[j].name): N_expected>0 but empty pool (cannot sample)")
    end
    esc = n_emit > 0 ? tot_drop / n_emit : 0.0
    esc <= 0.05 ||
        error(@sprintf("escaped fraction %.2f%% > 5%% — check the phantom/emitter co-registration",
                       100 * esc))
    s.t_meas_s > 0.0 || error("t_meas_s must be > 0")

    println("\n  OK: all sanity gates passed.")
end

main()
