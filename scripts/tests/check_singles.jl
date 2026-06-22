#!/usr/bin/env julia
# Validate a singles stack (scripts/simulate_source_mt.jl output), CSV or HDF5. Streams in
# bounded memory (CSV line-by-line; HDF5 in row-slices), asserts the STRUCTURAL invariants and
# reports the physics DISTRIBUTIONS, then exits nonzero on any hard-invariant violation so a
# run can be gated.
#
# Hard invariants (failure ⇒ exit 1):
#   • event_number non-decreasing (the streaming build_coincidences reader requires it);
#   • gamma ∈ {1,2}; iz ∈ [0,n_z); iphi ∈ [0,n_phi); nblocks ≥ 1; e_keV > 0.
# Reported: rows, events, reach %, gamma balance, phantom-scatter %, overspill %, energy + a
# coarse spectrum. The format is auto-detected from the file extension (.csv / .h5).
#
# TOML-config driven (reads prod/<tag>/singles.{csv,h5}, picking by [output].format):
#   julia --project=. scripts/check_singles.jl --config runs/sphere_water_csi.toml
#   julia --project=. scripts/check_singles.jl --config runs/sphere_water_csi.toml --singles path

using PTCryspMC
using ArgParse

const EBIN = 50.0          # energy histogram bin width [keV]
const NEBIN = 14           # 0 … 700 keV; the last bin also catches overflow
const E_HI = 600.0         # energy upper bound [keV] (511 + margin); 0 is valid (a near-zero
                           # forward-Compton recoil quantizes to 0 in 0.1 keV units)

function parse_cli()
    s = ArgParseSettings(description="Validate a singles stack (CSV or HDF5): structure + distributions.")
    @add_arg_table! s begin
        "--config";  help = "run config TOML"; required = true
        "--singles"; help = "override the singles path (default prod/<tag>/singles.{csv,h5})"; default = ""
    end
    parse_args(s)
end

# Running validation + distribution accumulator, fed one row at a time by either front-end.
mutable struct Stats
    rows::Int; maxev::Int; prev::Int
    g1::Int; g2::Int; scat::Int; over::Int
    esum::Float64; emin::Float64; emax::Float64
    ehist::Vector{Int}
    tsum::Float64; tmax::Float64
    v_order::Int; v_gamma::Int; v_block::Int; v_nblk::Int; v_energy::Int; v_time::Int
end
Stats() = Stats(0, 0, 0, 0, 0, 0, 0, 0.0, Inf, -Inf, zeros(Int, NEBIN), 0.0, -Inf, 0, 0, 0, 0, 0, 0)

@inline function feed!(st::Stats, ev::Int, g::Int, iz::Int, iphi::Int, nb::Int,
                       e::Float64, ph::Bool, t_rel::Float64, n_phi::Int, n_z::Int)
    st.rows += 1
    ev < st.prev                          && (st.v_order  += 1)
    (g == 1 || g == 2)                    || (st.v_gamma  += 1)
    (0 <= iz < n_z && 0 <= iphi < n_phi)  || (st.v_block  += 1)
    nb >= 1                               || (st.v_nblk   += 1)
    (0.0 <= e <= E_HI)                    || (st.v_energy += 1)
    (isfinite(t_rel) && t_rel >= 0.0)     || (st.v_time   += 1)   # t_rel = TOF + jitter ≥ 0
    st.tsum += t_rel; t_rel > st.tmax && (st.tmax = t_rel)
    st.prev = ev; ev > st.maxev && (st.maxev = ev)
    g == 1 ? (st.g1 += 1) : g == 2 && (st.g2 += 1)
    ph && (st.scat += 1)
    nb > 1 && (st.over += 1)
    st.esum += e; e < st.emin && (st.emin = e); e > st.emax && (st.emax = e)
    st.ehist[clamp(floor(Int, e / EBIN) + 1, 1, NEBIN)] += 1
end

function feed_csv!(st, path, n_phi, n_z)
    open(path, "r") do io
        header = split(strip(readline(io)), ',')
        col = Dict(String(h) => i for (i, h) in enumerate(header))
        for c in ("event_number", "gamma", "e_keV", "iz", "iphi", "nblocks", "phantom_scatter", "t_rel_ns")
            haskey(col, c) || error("singles stack is missing column '$c'")
        end
        iev = col["event_number"]; ig = col["gamma"]; ie = col["e_keV"]
        iiz = col["iz"]; iip = col["iphi"]; inb = col["nblocks"]; iph = col["phantom_scatter"]; it = col["t_rel_ns"]
        for line in eachline(io)
            isempty(line) && continue
            f = split(line, ',')
            feed!(st, parse(Int, f[iev]), parse(Int, f[ig]), parse(Int, f[iiz]), parse(Int, f[iip]),
                  parse(Int, f[inb]), parse(Float64, f[ie]), f[iph] == "1", parse(Float64, f[it]), n_phi, n_z)
        end
    end
end

function feed_hdf5!(st, path, n_phi, n_z)
    foreach_singles_hdf5(path) do b
        for i in 1:length(b)
            feed!(st, Int(b.event[i]), Int(b.gamma[i]), Int(b.iz[i]), Int(b.iphi[i]),
                  Int(b.nblocks[i]), decode_e(b.e[i]), b.phantom_scatter[i] == 1, Float64(b.t_rel[i]), n_phi, n_z)
        end
    end
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(prod_base(cfg)), tag)
    fmt = String(cfg_get(cfg, "output", "format", "csv"))
    default = joinpath(outdir, fmt == "hdf5" ? "singles.h5" : "singles.csv")
    singles = isempty(a["singles"]) ? default : rp(a["singles"])
    isfile(singles) || error("singles file '$singles' not found (run simulate_source_mt.jl first)")
    ishdf5 = endswith(singles, ".h5")

    # Scanner segmentation, to bound the block indices.
    datadir  = rp(cfg_get(cfg, "paths", "data", "data"))
    geomfile = rp(cfg_get(cfg, "geometry", "file", "geometry/geometry.json"))
    geom = load_geometry(geomfile, load_materials(datadir))
    geom.scanner === nothing && error("geometry $geomfile has no scanner section")
    n_phi = geom.scanner.n_phi; n_z = geom.scanner.n_z

    st = Stats()
    ishdf5 ? feed_hdf5!(st, singles, n_phi, n_z) : feed_csv!(st, singles, n_phi, n_z)

    nbad = st.v_order + st.v_gamma + st.v_block + st.v_nblk + st.v_energy + st.v_time
    ngamma = 2 * st.maxev      # max event ≈ nevents → total photons (a miss writes no row)

    println("singles: $singles  ($(ishdf5 ? "hdf5" : "csv"))")
    println("  rows (detected photons): $(st.rows)   events (max): $(st.maxev)   scanner $(n_phi)φ × $(n_z)z")
    if st.rows > 0
        println("  reach (rows / 2·events):   $(round(100*st.rows/max(ngamma,1), digits=1))%")
        println("  γ1 / γ2 balance:           $(st.g1) / $(st.g2)")
        println("  phantom-scattered:         $(round(100*st.scat/st.rows, digits=1))%")
        println("  overspill (nblocks>1):     $(round(100*st.over/st.rows, digits=1))%")
        println("  energy:  mean $(round(st.esum/st.rows, digits=1)) keV   range [$(round(st.emin,sigdigits=2)), $(round(st.emax,digits=1))] keV  (min ≈ 0: near-zero forward-Compton recoils, quantized)")
        println("  t_rel (TOF+jitter):        mean $(round(st.tsum/st.rows, digits=3)) ns   max $(round(st.tmax, digits=1)) ns")
        print("  spectrum (keV → %): ")
        for b in 1:NEBIN
            pct = round(100*st.ehist[b]/st.rows, digits=1)
            pct >= 0.1 && print("[$(Int((b-1)*EBIN))]$(pct)  ")
        end
        println()
    end

    println("invariants: order=$(st.v_order) gamma=$(st.v_gamma) block=$(st.v_block) nblocks=$(st.v_nblk) energy=$(st.v_energy) time=$(st.v_time)")
    if nbad == 0
        println("PASS ✓  ($(st.rows) rows, all invariants hold)")
    else
        println("FAIL ✗  ($nbad invariant violations)")
        exit(1)
    end
end

main()
