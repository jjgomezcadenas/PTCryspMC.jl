#!/usr/bin/env julia
# Build the RANDOM (accidental) LORs from a SINGLES stack — the third LOR category. Two singles
# from DIFFERENT annihilations whose ABSOLUTE arrival times fall within the coincidence window τ
# make a false coincidence. The singles already carry the per-photon time relative to their decay
# (`t_rel`, stamped by simulate_source_mt); this restores the absolute clock per single
# (`t_abs = event_time(event)·1e9 + t_rel`, in memory only), sorts, and pairs cross-event singles
# within τ (src/randoms.jl).
#
# Output: prod/<tag>/randoms.h5 — the LOR schema with truth=2, referenced to the EARLIER single's
# decay so the stored times stay small (t1 = t_rel[i], t2 = t_rel[i] + Δ, dT = Δ ≤ τ); dt_ns = NaN
# (no common emission point → no TOF correction). Detector smearing and the energy cut are NOT
# applied here — they come uniformly in the reco merge (lors_det.h5). So these are pre-energy-cut
# random candidates. The window τ is read from [timing].tau_ns (chosen from the DT study).
#
# Trues and randoms are NOT mutually exclusive: a single can appear in both a random here and its
# same-annihilation true (distinct LORs). This is the multiples approximation (no multiple-rejection,
# no opposition filter) — negligible (~randoms/trues), and the count matches 2τS². See pair_randoms.
#
#   julia --project=. scripts/build_randoms_from_singles.jl --config runs/sphere_water_csi.toml

using PTCryspMC
using ArgParse
using Printf

function parse_cli()
    s = ArgParseSettings(description="Build the random (accidental) LORs from a singles stack.")
    @add_arg_table! s begin
        "--config";  help = "run config TOML"; required = true
        "--singles"; help = "override the input singles path (default prod/<tag>/singles.{h5,csv})"; default = ""
        "--out";     help = "override the output randoms path (default prod/<tag>/randoms.h5)"; default = ""
    end
    parse_args(s)
end

# Columnar store of the CONTAINED singles (nblocks==1) kept in memory for the sort/pair. Positions
# stay quantized (Int16, as in the file) — decoded only when a pair is written.
struct RandSingles
    ev::Vector{Int32}
    t_abs::Vector{Float64}; t_rel::Vector{Float32}
    x::Vector{Int16}; y::Vector{Int16}; z::Vector{Int16}; e::Vector{Int16}
    iz::Vector{Int16}; iphi::Vector{Int16}; nscat::Vector{Int8}
    x0::Vector{Int16}; y0::Vector{Int16}; z0::Vector{Int16}
end
RandSingles() = RandSingles(Int32[], Float64[], Float32[], Int16[], Int16[], Int16[], Int16[],
                            Int16[], Int16[], Int8[], Int16[], Int16[], Int16[])
Base.length(r::RandSingles) = length(r.ev)

# `acts` is a per-isotope model vector (API) or length-1 (clinic/count, isotope always 0). The single's
# isotope column picks the model → its λ times the absolute-clock restore.
@inline function push_contained!(r::RandSingles, acts, ev, iso, x, y, z, e, iz, iphi, nb, x0, y0, z0, t_rel, nscat)
    nb == 1 || return                                    # contained hits only (same as the trues)
    t_abs = event_time(acts[iso + 1], ev) * 1.0e9 + t_rel # absolute clock [ns] with the isotope's λ
    push!(r.ev, Int32(ev)); push!(r.t_abs, t_abs); push!(r.t_rel, Float32(t_rel))
    push!(r.x, x); push!(r.y, y); push!(r.z, z); push!(r.e, e); push!(r.iz, iz); push!(r.iphi, iphi)
    push!(r.nscat, Int8(min(nscat, 127)))
    push!(r.x0, x0); push!(r.y0, y0); push!(r.z0, z0)
end

function collect_hdf5!(r::RandSingles, acts, path)
    foreach_singles_hdf5(path) do b
        for i in 1:length(b)
            push_contained!(r, acts, Int(b.event[i]), Int(b.isotope[i]), b.x[i], b.y[i], b.z[i], b.e[i],
                            b.iz[i], b.iphi[i], Int(b.nblocks[i]), b.x0[i], b.y0[i], b.z0[i],
                            Float64(b.t_rel[i]), Int(b.n_scatter[i]))
        end
    end
end

function collect_csv!(r::RandSingles, acts, path)
    open(path, "r") do io
        header = split(strip(readline(io)), ','); col = Dict(String(h) => i for (i, h) in enumerate(header))
        for c in ("event_number","gamma","x_mm","y_mm","z_mm","e_keV","iz","iphi","nblocks","n_scatter","isotope","x0_mm","y0_mm","z0_mm","t_rel_ns")
            haskey(col, c) || error("singles stack is missing column '$c'")
        end
        ie=col["event_number"]; ix=col["x_mm"]; iy=col["y_mm"]; iz=col["z_mm"]; iee=col["e_keV"]
        iiz=col["iz"]; iip=col["iphi"]; inb=col["nblocks"]; ins=col["n_scatter"]; iis=col["isotope"]
        ix0=col["x0_mm"]; iy0=col["y0_mm"]; iz0=col["z0_mm"]; it=col["t_rel_ns"]
        for line in eachline(io)
            isempty(line) && continue
            f = split(line, ',')
            push_contained!(r, acts, parse(Int, f[ie]), parse(Int, f[iis]),
                encode_xyz_mm(parse(Float64, f[ix])), encode_xyz_mm(parse(Float64, f[iy])), encode_xyz_mm(parse(Float64, f[iz])),
                encode_e_keV(parse(Float64, f[iee])), Int16(parse(Int, f[iiz])), Int16(parse(Int, f[iip])), parse(Int, f[inb]),
                encode_xyz_mm(parse(Float64, f[ix0])), encode_xyz_mm(parse(Float64, f[iy0])), encode_xyz_mm(parse(Float64, f[iz0])),
                parse(Float64, f[it]), parse(Int, f[ins]))
        end
    end
end

# Analytic random estimate: for a single rate S(t), expected cross-event pairs within τ ≈
# τ·∫S²dt ≈ τ·Σ n_k²/Δt over time bins (the uniform-Poisson forward-pair count over all contained
# singles). The same-event self-pairs it includes enter only at the negligible uniform within-τ
# weight, not the actual ns-clustered weight the measured pass excludes — so it tracks the measured
# cross-event randoms (ratio ≈ 1), but as an estimate, not an identity.
function analytic_randoms(t_abs::Vector{Float64}, tau::Float64; nbins::Int=200)
    lo, hi = extrema(t_abs); span = hi - lo
    span > 0 || return 0.0
    Δt = span / nbins; counts = zeros(Int, nbins)
    for t in t_abs
        counts[clamp(floor(Int, (t - lo) / Δt) + 1, 1, nbins)] += 1
    end
    tau * sum(c -> c^2 / Δt, counts)
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(prod_base(cfg)), tag)

    tau = Float64(cfg_get(cfg, "timing", "tau_ns", 0.0))
    tau > 0.0 || error("[timing].tau_ns must be set (> 0) — choose it from the DT study (examine_dt.jl)")
    crystal = String(cfg_get(cfg, "transport", "crystal_material", "CsI"))

    fmt = String(cfg_get(cfg, "output", "format", "csv"))
    singles = isempty(a["singles"]) ? joinpath(outdir, fmt == "hdf5" ? "singles.h5" : "singles.csv") : rp(a["singles"])
    isfile(singles) || error("singles file '$singles' not found (run simulate_source_mt.jl first)")
    ishdf5 = endswith(singles, ".h5")
    out = isempty(a["out"]) ? joinpath(outdir, "randoms.h5") : rp(a["out"])

    # Timing: per-isotope when the singles carry `isotope_half_lives` (API — one truncated-exponential
    # model per isotope over [0, t_meas]), else the single config `ActivityModel` (clinic/count; the
    # isotope column is all 0, so acts[1] is used for every single).
    hls = ishdf5 ? Vector{Float64}(singles_hdf5_attr(singles, "isotope_half_lives", Float64[])) : Float64[]
    if !isempty(hls)
        tmeas = Float64(singles_hdf5_attr(singles, "t_meas_s", 600.0))
        tseed = Int(singles_hdf5_attr(singles, "time_seed", 1234))
        acts  = ActivityModel[ActivityModel(; t0=0.0, t1=tmeas, half_life_s=h, seed=tseed) for h in hls]
        tdesc = "per-isotope ($(length(acts)) isotopes, t_meas=$(tmeas)s)"
    else
        acts  = ActivityModel[ActivityModel(cfg)]
        tdesc = "single-isotope (t½=$(round(log(2.0)/acts[1].λ, digits=1))s, window [$(acts[1].t0),$(acts[1].t1)]s)"
    end

    println("run '$tag' randoms: τ=$tau ns, crystal $crystal, timing $tdesc")
    println("  singles ($(ishdf5 ? "hdf5" : "csv")): $singles")

    r = RandSingles()
    ishdf5 ? collect_hdf5!(r, acts, singles) : collect_csv!(r, acts, singles)
    println("  contained singles kept: $(length(r))")

    meta = Dict{String,Any}("scenario_tag" => tag, "mode" => "randoms", "has_randoms" => true,
        "crystal" => crystal, "tau_ns" => tau, "t_relative_to_decay" => true,
        "t0_s" => acts[1].t0, "t1_s" => acts[1].t1, "n_isotopes" => length(acts), "time_seed" => Int(acts[1].seed))
    w = CoincidenceWriter(out, meta)

    # Emit one random LOR per cross-event pair (i earlier, j later), referenced to i's decay.
    emit = function (i, j, Δ)
        push_coincidence!(w, r.ev[i],
            decode_xyz(r.x[i]), decode_xyz(r.y[i]), decode_xyz(r.z[i]), decode_e(r.e[i]), r.t_rel[i], r.iz[i], r.iphi[i], Int(r.nscat[i]),
            decode_xyz(r.x[j]), decode_xyz(r.y[j]), decode_xyz(r.z[j]), decode_e(r.e[j]), r.t_rel[i] + Δ, r.iz[j], r.iphi[j], Int(r.nscat[j]),
            NaN32, decode_xyz(r.x0[i]), decode_xyz(r.y0[i]), decode_xyz(r.z0[i]), TRUTH_RANDOM)
    end
    nrand = pair_randoms(emit, r.t_abs, r.ev, tau)

    nevents = ishdf5 ? Int(singles_hdf5_attr(singles, "nevents", 0)) : 0
    set_lor_attr!(w, "nevents", nevents)
    close(w)

    expect = analytic_randoms(r.t_abs, tau)
    println("wrote $nrand random LORs -> $out")
    @printf "  analytic τ·∫S²dt ≈ %.0f   (measured/analytic = %.2f)\n" expect nrand/max(expect, 1e-9)
    println("  (truth=2; pre-energy-cut candidates — smearing/energy/DT cut applied in the reco merge)")
end

main()
