#!/usr/bin/env julia
# Validate a LOR (coincidence) HDF5 file (build_coincidences_from_singles.jl output). Streams it
# in row-slices, asserts the STRUCTURAL invariants and reports the physics DISTRIBUTIONS, then
# exits nonzero on any hard-invariant violation so a run can be gated.
#
# Hard invariants (failure ⇒ exit 1):
#   • event non-decreasing; truth ∈ {0,1,2}; both hits' iz ∈ [0,n_z), iphi ∈ [0,n_phi);
#   • 0 ≤ e1,e2 ≤ 600 keV.
# Reported: rows, acceptance %, truth split (true/scatter/random), hit energies, and the block
# OPPOSITION (a coincidence should join two roughly opposite crystals — Δφ ≈ n_phi/2).
#
# TOML-config driven (reads prod/<tag>/lors_{truth|det}.h5, mode from [detector]):
#   julia --project=. scripts/tests/check_lors.jl --config runs/sphere_water_csi.toml
#   julia --project=. scripts/tests/check_lors.jl --config runs/sphere_water_csi.toml --lors path.h5

using PTCryspMC
using ArgParse

const E_HI = 600.0         # energy upper bound [keV] (511 + margin)

function parse_cli()
    s = ArgParseSettings(description="Validate a LOR HDF5 file: structure + distributions.")
    @add_arg_table! s begin
        "--config"; help = "run config TOML"; required = true
        "--lors";   help = "override the LOR path (default prod/<tag>/lors_{truth|det}.h5)"; default = ""
    end
    parse_args(s)
end

mutable struct Stats
    rows::Int; prev::Int; maxev::Int
    ntrue::Int; nscat::Int; nrand::Int
    e1sum::Float64; e2sum::Float64
    opp::Int                                  # hits roughly opposite in φ
    v_order::Int; v_truth::Int; v_block::Int; v_energy::Int
end
Stats() = Stats(0, 0, 0, 0, 0, 0, 0.0, 0.0, 0, 0, 0, 0, 0)

@inline function feed!(st::Stats, ev::Int, truth::Int, iz1::Int, iphi1::Int, iz2::Int, iphi2::Int,
                       e1::Float64, e2::Float64, n_phi::Int, n_z::Int)
    st.rows += 1
    ev < st.prev                                          && (st.v_order  += 1)
    (0 <= truth <= 2)                                     || (st.v_truth  += 1)
    (0 <= iz1 < n_z && 0 <= iphi1 < n_phi &&
     0 <= iz2 < n_z && 0 <= iphi2 < n_phi)                || (st.v_block  += 1)
    (0.0 <= e1 <= E_HI && 0.0 <= e2 <= E_HI)              || (st.v_energy += 1)
    st.prev = ev; ev > st.maxev && (st.maxev = ev)
    truth == 0 ? (st.ntrue += 1) : truth == 1 ? (st.nscat += 1) : (st.nrand += 1)
    st.e1sum += e1; st.e2sum += e2
    dφ = abs(iphi1 - iphi2); circ = min(dφ, n_phi - dφ)   # circular block separation
    abs(circ - n_phi ÷ 2) <= n_phi ÷ 8 && (st.opp += 1)   # within ±¼ of opposite
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(prod_base(cfg)), tag)
    sigma_xyz = Float64(cfg_get(cfg, "detector", "sigma_xyz_mm", 0.0))
    eres   = Float64(cfg_get(cfg, "detector", "eres", 0.0))
    emin   = Float64(cfg_get(cfg, "detector", "emin_keV", 0.0))
    window = Float64(cfg_get(cfg, "detector", "window_fwhm", 0.0))
    mode = (sigma_xyz > 0.0 || eres > 0.0 || window > 0.0 || emin > 0.0) ? "det" : "truth"
    lors = isempty(a["lors"]) ? joinpath(outdir, "lors_$mode.h5") : rp(a["lors"])
    isfile(lors) || error("LOR file '$lors' not found (run build_coincidences_from_singles.jl first)")

    datadir  = rp(cfg_get(cfg, "paths", "data", "data"))
    geomfile = rp(cfg_get(cfg, "geometry", "file", "geometry/geometry.json"))
    geom = load_geometry(geomfile, load_materials(datadir))
    geom.scanner === nothing && error("geometry $geomfile has no scanner section")
    n_phi = geom.scanner.n_phi; n_z = geom.scanner.n_z

    nevents     = Int(singles_hdf5_attr(lors, "nevents", 0))
    has_randoms = Bool(singles_hdf5_attr(lors, "has_randoms", false))

    st = Stats()
    foreach_coincidences_hdf5(lors) do b
        for i in 1:length(b)
            feed!(st, Int(b.event[i]), Int(b.truth[i]), Int(b.iz1[i]), Int(b.iphi1[i]),
                  Int(b.iz2[i]), Int(b.iphi2[i]), decode_e(b.e1[i]), decode_e(b.e2[i]), n_phi, n_z)
        end
    end

    nbad = st.v_order + st.v_truth + st.v_block + st.v_energy
    println("lors: $lors  (mode $mode, has_randoms=$has_randoms)")
    println("  LORs: $(st.rows)   events: $nevents   scanner $(n_phi)φ × $(n_z)z")
    if st.rows > 0
        nevents > 0 && println("  acceptance:    $(round(100*st.rows/nevents, digits=1))% of annihilations")
        println("  truth split:   $(st.ntrue) true / $(st.nscat) scatter / $(st.nrand) random " *
                "($(round(100*st.ntrue/st.rows,digits=1))% / $(round(100*st.nscat/st.rows,digits=1))% / $(round(100*st.nrand/st.rows,digits=1))%)")
        println("  hit energy:    mean e1 $(round(st.e1sum/st.rows,digits=1)) keV, e2 $(round(st.e2sum/st.rows,digits=1)) keV")
        println("  opposition:    $(round(100*st.opp/st.rows,digits=1))% of LORs join roughly opposite blocks (Δφ ≈ n_phi/2)")
    end
    println("invariants: order=$(st.v_order) truth=$(st.v_truth) block=$(st.v_block) energy=$(st.v_energy)")
    if nbad == 0
        println("PASS ✓  ($(st.rows) LORs, all invariants hold)")
    else
        println("FAIL ✗  ($nbad invariant violations)")
        exit(1)
    end
end

main()
