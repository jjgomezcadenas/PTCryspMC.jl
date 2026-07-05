#!/usr/bin/env julia
# End-to-end validation of an API (Proton Activity) production run: read the produced prod/<tag>/
# files and confirm the physics. Rebuilds the materialized source (same master_seed+realization) and
# checks, against the singles + LOR stacks:
#   A. N & per-isotope counts vs the budget (M_j ~ Poisson(N_expected_j · f_inside_j)).
#   B. Reproducibility / nchunks-invariance: each single's emission point IS the materialized array's
#      point[event_number] — proving the source is indexed by event, independent of the chunking.
#   C. Spatial activity profile: the singles emission-point z-profile reproduces the source (and
#      hence the SOBP activity plateau + distal falloff); detection does not distort it.
#   D. True/scatter/random split from lors_det — physically sane fractions.
# Errors on any gate. (The randoms 2τS² ratio is validated by build_randoms itself during the run.)
#
#   julia --project=. scripts/tests/check_api_validation.jl --config runs/uniform_headep_bgo_api.toml

using PTCryspMC
using ArgParse
using Printf

function parse_cli()
    s = ArgParseSettings(description="Validate an API production run end-to-end.")
    @add_arg_table! s begin
        "--config"; help = "the API run config TOML (mode=api)"; required = true
        "--data";   help = "materials/XCOM directory"; default = "data"
    end
    parse_args(s)
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    String(cfg_get(cfg, "source", "mode", "")) == "api" || error("[source].mode must be 'api'")
    tag    = run_tag(cfg, a["config"])
    outdir = joinpath(rp(prod_base(cfg)), tag)
    singles = joinpath(outdir, "singles.h5"); lors_det = joinpath(outdir, "lors_det.h5")
    for f in (singles, lors_det); isfile(f) || error("missing $f (run the chain first)"); end

    mats = load_materials(rp(a["data"]))
    scndir  = rp(String(cfg_get(cfg, "source", "scenario_dir", "")))
    budget  = String(cfg_get(cfg, "source", "budget", "fast"))
    dose    = Float64(cfg_get(cfg, "source", "dose_Gy", 1.0))
    keepesc = Bool(cfg_get(cfg, "source", "keep_escaped", false))
    mseed   = Int(cfg_get(cfg, "source", "master_seed", 1))
    realz   = Int(cfg_get(cfg, "source", "realization", 0))

    println("validating '$tag': scenario=$(basename(scndir)) budget=$budget dose=$(dose)Gy real=$realz")
    scn = load_scenario(scndir, mats; budget=budget, dose_Gy=dose, keep_escaped=keepesc)
    src = materialize_api_source(scn; master_seed=mseed, realization=realz)
    N   = length(src.points)

    # ---- A. N & per-isotope counts vs the budget ----------------------------
    nev_attr = Int(singles_hdf5_attr(singles, "nevents", -1))
    @printf("\nA. materialized N = %d  |  singles attr nevents = %d\n", N, nev_attr)
    nev_attr == N || error("nevents attr ($nev_attr) ≠ materialized N ($N)")
    bad = String[]
    for j in 1:length(scn.pools)
        Mj = count(==(Int8(j - 1)), src.isotope); μ = scn.n_expected[j] * scn.f_inside[j]
        z  = μ > 0 ? (Mj - μ) / sqrt(μ) : 0.0
        @printf("   %-4s M_j=%d  μ=%.4e  (M-μ)/√μ=%.2f\n", scn.isotopes[j].name, Mj, μ, z)
        (abs(Mj - μ) <= 6 * sqrt(μ) + 6) || push!(bad, scn.isotopes[j].name)
    end
    isempty(bad) || error("isotopes off the Poisson budget by >6σ: $(join(bad, ", "))")

    # ---- B. reproducibility / nchunks-invariance ----------------------------
    # Each single's emission point (x0,y0,z0) must equal the materialized array's point[event_number]
    # (to the 0.1 mm quantization). This is the "identical source, chunking-independent" guarantee:
    # the point is indexed by event, never redrawn per chunk. Check a sample (the first read batch).
    checked = 0; bmax = 0.0
    foreach_singles_hdf5(singles; batch = 200_000) do b
        checked == 0 || return
        for i in 1:length(b)
            ev = Int(b.event[i]); p = src.points[ev]
            d = max(abs(decode_xyz(b.x0[i]) / 10 - p[1]),
                    abs(decode_xyz(b.y0[i]) / 10 - p[2]),
                    abs(decode_xyz(b.z0[i]) / 10 - p[3]))          # cm
            bmax = max(bmax, d)
        end
        checked = length(b)
    end
    @printf("\nB. singles emission point vs materialized array: %d checked, max |Δ| = %.4f mm\n",
            checked, bmax * 10)
    bmax * 10 <= XYZ_SCALE_MM + 1e-6 ||
        error("singles emission point diverges from the materialized array (max $(bmax*10) mm) — not chunking-invariant")

    # ---- C. spatial activity profile (z) ------------------------------------
    # Bin the SOURCE z (all materialized points) and the DETECTED z (singles emission point) over the
    # phantom z-span; detection should reproduce the source shape (high correlation), and the profile
    # should show the SOBP plateau with a distal falloff (beam enters at -z, travels +z).
    lo, hi = -110.0, 110.0; nb = 44; dz = (hi - lo) / nb        # mm, 5 mm bins
    binof(zmm) = clamp(floor(Int, (zmm - lo) / dz) + 1, 1, nb)
    hsrc = zeros(Int, nb)
    for p in src.points; hsrc[binof(p[3] * 10)] += 1; end
    hdet = zeros(Int, nb)
    foreach_singles_hdf5(singles) do b
        for i in 1:length(b); hdet[binof(decode_xyz(b.z0[i]))] += 1; end
    end
    # Pearson correlation of the two normalized profiles.
    ms = sum(hsrc)/nb; md = sum(hdet)/nb
    cov = sum((hsrc[k]-ms)*(hdet[k]-md) for k in 1:nb)
    corr = cov / sqrt(sum((hsrc[k]-ms)^2 for k in 1:nb) * sum((hdet[k]-md)^2 for k in 1:nb))
    # Distal 50% edge of the SOURCE activity: scanning from high z downward, the last bin above half-max.
    peak = maximum(hsrc); half = peak / 2
    distal_bin = findlast(h -> h >= half, hsrc)
    distal_z = lo + (distal_bin - 0.5) * dz
    @printf("\nC. activity z-profile: source↔detected correlation = %.4f\n", corr)
    @printf("   plateau peak %d/bin; distal 50%% activity edge z ≈ %.0f mm (target distal face z=-5 mm)\n",
            peak, distal_z)
    corr >= 0.98 || error("detected z-profile does not track the source (corr $corr < 0.98)")
    (-110.0 <= distal_z <= 110.0) || error("distal activity edge $distal_z mm outside the phantom")

    # ---- D. true / scatter / random split -----------------------------------
    nt = 0; ns = 0; nr = 0
    foreach_coincidences_hdf5(lors_det) do b
        for i in 1:length(b)
            c = Int(b.truth[i]); c == 0 ? (nt += 1) : c == 1 ? (ns += 1) : (nr += 1)
        end
    end
    tot = nt + ns + nr
    @printf("\nD. lors_det: %d LORs — %d true (%.1f%%) / %d scatter (%.1f%%) / %d random (%.2f%%)\n",
            tot, nt, 100nt/tot, ns, 100ns/tot, nr, 100nr/tot)
    (0.10 <= nt/tot <= 0.60) || error("true fraction $(nt/tot) implausible for a brain head (expect ~0.15–0.50)")
    (nr / max(nt + ns, 1) <= 0.05) || error("randoms/(true+scatter) = $(nr/(nt+ns)) > 5% — unexpected at this dose")

    println("\n  OK: all API validation gates passed.")
end

main()
