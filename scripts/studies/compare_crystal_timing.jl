#!/usr/bin/env julia
# Compare the coincidence timing across crystals — why the coincidence-window variable dT =
# |t1−t2| comes out similar for CsI and BGO even though BGO's lower light yield gives a much
# worse single-photon jitter. For each run config it reports:
#
#   • the ANALYTIC single-photon jitter from the crystal DB: 1/(N_det·r0), N_det = yield·E·pde,
#     r0 = Σ wₖ/τₖ — the per-deposit timing fluctuation (∝ 1/√… no: ∝ 1/N_det), at a few energies;
#   • the MEASURED dT distribution from lors_truth.h5 (p50/p99/p99.9 — the window lives in the
#     tail), the residual dt = |Δt0|−TOF_diff core (the intrinsic timing resolution, TOF removed),
#     and the scatter fraction (a photofraction proxy).
#
# The point: dT's TAIL (hence the window τ) is similar because BGO's lower yield (worse) is offset
# by its higher photofraction (more full-energy deposits → fewer bad-timing low-energy ones), while
# the CORE / residual still shows BGO's genuinely worse timing.
#
#   julia --project=. scripts/studies/compare_crystal_timing.jl \
#       --configs runs/sphere_water_csi.toml runs/sphere_water_bgo.toml

using PTCryspMC
using ArgParse
using Printf

function parse_cli()
    s = ArgParseSettings(description="Compare coincidence timing (dT, residual, jitter) across crystals.")
    @add_arg_table! s begin
        "--configs"; help = "one or more run config TOMLs"; nargs = '+'; required = true
    end
    parse_args(s)
end

# Mean single-photon jitter [ns] = 1/(N_det·r0): the analytic mean of first_photon_jitter.
function jitter_mean(mat, E_MeV)
    N = mat.light_yield * E_MeV * mat.pde
    r0 = 0.0
    for k in eachindex(mat.scint_decay_ns)
        mat.scint_decay_ns[k] > 0 && (r0 += mat.scint_decay_w[k] / mat.scint_decay_ns[k])
    end
    (N > 0 && r0 > 0) ? 1.0 / (N * r0) : 0.0
end

pct(v, q) = (n = length(v); n == 0 ? NaN : v[clamp(ceil(Int, q * n), 1, n)])  # v must be sorted
mean(v) = isempty(v) ? NaN : sum(v) / length(v)
stdev(v) = (n = length(v); n < 2 ? NaN : (m = mean(v); sqrt(sum(x -> (x - m)^2, v) / (n - 1))))

struct CrystalTiming
    tag::String; crystal::String
    yield::Float64; decay::Vector{Float64}; pde::Float64; jit511::Float64
    nlor::Int; scatfrac::Float64
    dT50::Float64; dT99::Float64; dT999::Float64
    res_med::Float64; res_core_std::Float64
end

function analyse(cfg_path, REPO)
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))
    cfg = read_config(cfg_path)
    tag = run_tag(cfg, cfg_path)
    crystal = String(cfg_get(cfg, "transport", "crystal_material", "CsI"))
    mat = load_material(joinpath(REPO, "data"), crystal)
    lors = joinpath(rp(prod_base(cfg)), tag, "lors_truth.h5")
    isfile(lors) || error("$lors not found (run build_true_coincidences_from_singles.jl first)")

    dTs = Float32[]; absdt = Float32[]; nscat = 0
    foreach_coincidences_hdf5(lors) do b
        for i in 1:length(b)
            push!(dTs, abs(b.t1[i] - b.t2[i]))
            push!(absdt, abs(b.dt[i]))
            b.truth[i] == 1 && (nscat += 1)
        end
    end
    sort!(dTs)
    core = filter(<(2.0f0), absdt)                         # exclude the low-energy jitter tail
    CrystalTiming(tag, crystal, mat.light_yield, mat.scint_decay_ns, mat.pde,
                  jitter_mean(mat, 0.511), length(dTs), nscat / max(length(dTs), 1),
                  pct(dTs, 0.50), pct(dTs, 0.99), pct(dTs, 0.999),
                  pct(sort(absdt), 0.50), stdev(core))
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rows = [analyse(c, REPO) for c in a["configs"]]

    println("Single-photon jitter 1/(N_det·r0) [ns] vs deposited energy (from the crystal DB):")
    @printf "  %-6s %12s %16s   %8s %8s %8s %8s\n" "crystal" "yield[γ/MeV]" "decay[ns]" "511keV" "200keV" "100keV" "50keV"
    for r in rows
        mat = load_material(joinpath(REPO, "data"), r.crystal)
        dstr = join(Int.(round.(r.decay)), "/")
        @printf "  %-6s %12.2g %16s   %8.3f %8.3f %8.3f %8.3f\n" r.crystal r.yield dstr jitter_mean(mat,0.511) jitter_mean(mat,0.200) jitter_mean(mat,0.100) jitter_mean(mat,0.050)
    end

    println("\nMeasured coincidence timing from lors_truth.h5:")
    @printf "  %-18s %8s %8s | %8s %8s %8s | %10s %10s\n" "tag" "LORs" "scatter" "dT_p50" "dT_p99" "dT_p99.9" "res_med|dt|" "res_core_std"
    @printf "  %-18s %8s %8s | %8s %8s %8s | %10s %10s\n" "" "" "%" "[ns]" "[ns]" "[ns]" "[ns]" "[ns]"
    for r in rows
        nlab = @sprintf("%.2fM", r.nlor/1e6)
        @printf "  %-18s %8s %7.1f%% | %8.3f %8.3f %8.2f | %10.3f %10.3f\n" r.tag nlab 100*r.scatfrac r.dT50 r.dT99 r.dT999 r.res_med r.res_core_std
    end

    println("\nReading: the dT TAIL (p99/p99.9 → the window τ) is similar across crystals — lower yield")
    println("is offset by higher photofraction (less scatter). The residual core shows the intrinsic")
    println("timing resolution, where the lower-yield crystal is genuinely worse.")
end

main()
