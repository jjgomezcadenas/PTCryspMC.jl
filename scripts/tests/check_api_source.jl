#!/usr/bin/env julia
# QA for the API source phase-1 materialization (`materialize_api_source`): draw one realization of
# the annihilation-point array from a frozen scenario and check it against the budget — the per-
# isotope M_j ~ Poisson(N_expected_j · f_inside_j), the total ΣM_j, reproducibility per
# (master_seed, realization) (the "identical source" contract, independent of the transport
# chunking), that a different realization differs, and that every materialized point sits inside the
# phantom (the escaped positrons were dropped at read). Runs sanity gates and errors on violation.
#
#   julia --project=. scripts/tests/check_api_source.jl --scenario ~/Projects/ptcrysp-scenarios/scenarios/uniform_headep_sobp_1e8
#   julia --project=. scripts/tests/check_api_source.jl --scenario <dir> --budget offline --dose 2 --realization 3

using PTCryspMC
using ArgParse
using Printf

function parse_cli()
    s = ArgParseSettings(description="Check the API source phase-1 Poisson materialization.")
    @add_arg_table! s begin
        "--scenario";     help = "scenario directory"; required = true
        "--budget";       help = "timing budget: fast | inroom | offline"; default = "fast"
        "--dose";         help = "clinical dose [Gy]"; arg_type = Float64; default = 1.0
        "--master-seed";  help = "master seed for the source draw"; arg_type = Int; default = 1
        "--realization";  help = "realization index (with master_seed, fixes the source)"; arg_type = Int; default = 0
        "--keep-escaped"; help = "keep escaped positrons instead of dropping"; action = :store_true
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
    scn  = load_scenario(rp(a["scenario"]), mats;
                         budget=a["budget"], dose_Gy=a["dose"], keep_escaped=a["keep-escaped"])
    seed = a["master-seed"]; real = a["realization"]

    t = @elapsed src = materialize_api_source(scn; master_seed=seed, realization=real)
    N = length(src.points)
    @printf("scenario '%s'  budget=%s  dose=%.3g Gy  seed=%d  realization=%d\n",
            scn.name, scn.budget, scn.dose_Gy, seed, real)
    @printf("materialized N = %d events in %.1fs  (~%.2f GB for points)\n", N, t, N * 24 / 1e9)

    @printf("\n  %-5s %14s %14s %9s %9s\n", "iso", "M_j", "μ=N·f_in", "M/μ", "(M-μ)/√μ")
    bad = String[]
    tot_mu = 0.0
    for j in 1:length(scn.pools)
        Mj = count(==(Int8(j - 1)), src.isotope)
        μ  = scn.n_expected[j] * scn.f_inside[j]; tot_mu += μ
        z  = μ > 0 ? (Mj - μ) / sqrt(μ) : 0.0
        @printf("  %-5s %14d %14.4e %9.4f %9.2f\n", scn.isotopes[j].name, Mj, μ,
                μ > 0 ? Mj / μ : NaN, z)
        # 6σ + small slack: a genuine Poisson draw is essentially never this far out.
        (abs(Mj - μ) <= 6 * sqrt(μ) + 6) || push!(bad, scn.isotopes[j].name)
    end
    @printf("  %-5s %14d %14.4e %9.4f\n", "all", N, tot_mu, tot_mu > 0 ? N / tot_mu : NaN)

    # Reproducibility (the "identical source" contract) and realization independence.
    src2 = materialize_api_source(scn; master_seed=seed, realization=real)
    reproducible = (src2.points == src.points) && (src2.isotope == src.isotope)
    src3 = materialize_api_source(scn; master_seed=seed, realization=real + 1)
    differs = (src3.points != src.points)

    # Containment: a materialized point must be inside the phantom (escaped were dropped at read).
    nsamp = min(N, 2_000_000)
    inside = count(i -> is_inside(scn.phantom, src.points[i]), 1:nsamp)
    @printf("\n  reproducible (same seed+realization): %s\n", reproducible)
    @printf("  realization %d differs from %d:         %s\n", real + 1, real, differs)
    @printf("  points inside phantom (sample %d):   %d\n", nsamp, inside)

    # ---- sanity gates -------------------------------------------------------
    isempty(bad) || error("isotopes off the Poisson budget by >6σ: $(join(bad, ", "))")
    reproducible || error("materialization not reproducible for the same (master_seed, realization)")
    differs || error("realization $(real+1) is identical to $real — the source is not re-realizing")
    (a["keep-escaped"] || inside == nsamp) ||
        error("$(nsamp - inside) of $nsamp materialized points fall OUTSIDE the phantom (escaped not dropped?)")

    println("\n  OK: all sanity gates passed.")
end

main()
