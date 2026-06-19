#!/usr/bin/env julia
# Validate a singles stack (scripts/simulate_source_mt.jl output). Streams the CSV line by line
# in O(1) memory (the stack can be > 1 GB), asserts the STRUCTURAL invariants and reports the
# physics DISTRIBUTIONS, then exits nonzero on any hard-invariant violation so a run can be
# gated in a pipeline.
#
# Hard invariants (failure ⇒ exit 1):
#   • the expected header columns are present;
#   • event_number is non-decreasing (the streaming build_coincidences reader requires it);
#   • gamma ∈ {1,2}; iz ∈ [0,n_z); iphi ∈ [0,n_phi); nblocks ≥ 1; e_keV > 0.
# Reported (informational): rows, events, reach %, phantom-scatter %, overspill %, the gamma
# balance, mean energy and a coarse energy spectrum.
#
# TOML-config driven (reads output/<tag>/singles.csv), mirroring simulate_source_mt.jl:
#   julia --project=. scripts/check_singles.jl --config runs/sphere_water_csi.toml
#   julia --project=. scripts/check_singles.jl --config runs/sphere_water_csi.toml --singles path/to.csv

using PTCryspMC
using ArgParse

function parse_cli()
    s = ArgParseSettings(description="Validate a singles stack (structure + distributions).")
    @add_arg_table! s begin
        "--config";  help = "run config TOML"; required = true
        "--singles"; help = "override the singles path (default output/<tag>/singles.csv)"; default = ""
    end
    parse_args(s)
end

const EBIN = 50.0          # energy histogram bin width [keV]
const NEBIN = 14           # bins span 0 … 700 keV; the last bin also catches overflow

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(cfg_get(cfg, "output", "dir", "output")), tag)
    singles = isempty(a["singles"]) ? joinpath(outdir, "singles.csv") : rp(a["singles"])
    isfile(singles) || error("singles file '$singles' not found (run simulate_source_mt.jl first)")

    # Scanner segmentation, to bound the block indices.
    datadir  = rp(cfg_get(cfg, "paths", "data", "data"))
    geomfile = rp(cfg_get(cfg, "geometry", "file", "geometry/geometry.json"))
    geom = load_geometry(geomfile, load_materials(datadir))
    geom.scanner === nothing && error("geometry $geomfile has no scanner section")
    n_phi = geom.scanner.n_phi; n_z = geom.scanner.n_z

    # Counters (O(1) memory).
    rows = 0; maxev = 0; prev = 0
    g1 = 0; g2 = 0; scat = 0; over = 0
    esum = 0.0; emin = Inf; emax = -Inf
    ehist = zeros(Int, NEBIN)
    v_order = 0; v_gamma = 0; v_block = 0; v_nblk = 0; v_energy = 0

    open(singles, "r") do io
        header = split(strip(readline(io)), ',')
        col = Dict(String(h) => i for (i, h) in enumerate(header))
        for c in ("event_number", "gamma", "e_keV", "iz", "iphi", "nblocks", "phantom_scatter")
            haskey(col, c) || error("singles stack is missing column '$c'")
        end
        iev = col["event_number"]; ig = col["gamma"]; ie = col["e_keV"]
        iiz = col["iz"]; iip = col["iphi"]; inb = col["nblocks"]; iph = col["phantom_scatter"]

        for line in eachline(io)
            isempty(line) && continue
            f = split(line, ',')
            ev = parse(Int, f[iev]); g = parse(Int, f[ig])
            iz = parse(Int, f[iiz]); ip = parse(Int, f[iip]); nb = parse(Int, f[inb])
            e  = parse(Float64, f[ie])
            rows += 1

            ev < prev      && (v_order  += 1)
            (g == 1 || g == 2) || (v_gamma += 1)
            (0 <= iz < n_z && 0 <= ip < n_phi) || (v_block += 1)
            nb >= 1        || (v_nblk   += 1)
            e  > 0.0       || (v_energy += 1)
            prev = ev; ev > maxev && (maxev = ev)

            g == 1 ? (g1 += 1) : g == 2 && (g2 += 1)
            f[iph] == "1" && (scat += 1)
            nb > 1 && (over += 1)
            esum += e; e < emin && (emin = e); e > emax && (emax = e)
            b = clamp(floor(Int, e / EBIN) + 1, 1, NEBIN); ehist[b] += 1
        end
    end

    nbad = v_order + v_gamma + v_block + v_nblk + v_energy
    ngamma = 2 * maxev      # max event ≈ nevents → total photons (a miss writes no row)

    println("singles: $singles")
    println("  rows (detected photons): $rows   events (max): $maxev   scanner $(n_phi)φ × $(n_z)z")
    if rows > 0
        println("  reach (rows / 2·events):   $(round(100*rows/max(ngamma,1), digits=1))%")
        println("  γ1 / γ2 balance:           $g1 / $g2")
        println("  phantom-scattered:         $(round(100*scat/rows, digits=1))%")
        println("  overspill (nblocks>1):     $(round(100*over/rows, digits=1))%")
        println("  energy:  mean $(round(esum/rows, digits=1)) keV   range [$(round(emin,sigdigits=2)), $(round(emax,digits=1))] keV  (tiny min = near-zero forward-Compton recoils)")
        print("  spectrum (keV → %): ")
        for b in 1:NEBIN
            pct = round(100*ehist[b]/rows, digits=1)
            pct >= 0.1 && print("[$(Int((b-1)*EBIN))]$(pct)  ")
        end
        println()
    end

    println("invariants: order=$v_order gamma=$v_gamma block=$v_block nblocks=$v_nblk energy=$v_energy")
    if nbad == 0
        println("PASS ✓  ($rows rows, all invariants hold)")
    else
        println("FAIL ✗  ($nbad invariant violations)")
        exit(1)
    end
end

main()
