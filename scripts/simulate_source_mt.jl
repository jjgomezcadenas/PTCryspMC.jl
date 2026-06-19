#!/usr/bin/env julia
# Production singles driver — TOML-config driven, MULTI-THREADED, singles-only output.
#
# Same source + transport as scripts/simulate_phantom.jl (a back-to-back 511 keV source fills
# the phantom and emits pairs with acollinearity, each photon navigated phantom → air → ring),
# but it writes the REDUCED *singles* stack — one row per detected photon — via the allocation-
# free `navigate_single_photons`, so it scales to the ~10⁸-decay production runs. The full per-
# interaction stack (CNN R&D, detailed plots) stays in simulate_phantom.jl.
#
# Parallelism: the N events are split into `nchunks` CONTIGUOUS ranges (default 8·nthreads).
# Each chunk runs on a thread with its OWN pre-allocated RNG (`MersenneTwister(seed + chunk-1)`,
# chunk 1 = base seed) and streams its rows to its own part-file `singles.part<c>.csv`. After
# the parallel region the part-files are concatenated IN CHUNK ORDER into `singles.csv`, so the
# output is globally event-ordered (the streaming build_coincidences reader needs no change).
# Reproducible over (seed, nchunks, N), independent of thread count / scheduling.
#
# Singles schema (one row per photon that deposited in the ring; a miss writes nothing):
#   event_number, gamma, x_mm, y_mm, z_mm, e_keV, iz, iphi, nblocks, phantom_scatter, x0_mm, y0_mm, z0_mm
# x,y,z / iz,iphi = first crystal interaction (LOR point + block); e_keV = summed block energy
# (truth, unsmeared — smearing stays in build_coincidences); nblocks = distinct blocks (1 =
# contained, >1 = overspill). Times are not here yet (Step 5, randoms).
#
# Run from the repo root (note -t for threads):
#   julia -t 18 --project=. scripts/simulate_source_mt.jl --config runs/sphere_water_csi.toml
#   julia -t 1  --project=. scripts/simulate_source_mt.jl --config runs/sphere_water_csi.toml --nchunks 1

using PTCryspMC
using ArgParse
using Random
using Base.Threads: @threads, nthreads

function parse_cli()
    s = ArgParseSettings(description="Production singles driver (multi-threaded, TOML-config driven).")
    @add_arg_table! s begin
        "--config";  help = "run config TOML"; required = true
        "--nevents"; help = "override [source].nevents"; arg_type = Int; default = -1
        "--seed";    help = "override [transport].seed"; arg_type = Int; default = -1
        "--nchunks"; help = "number of event chunks (default 8·nthreads)"; arg_type = Int; default = -1
        "--format";  help = "override [output].format (csv | hdf5)"; default = ""
    end
    parse_args(s)
end

# Write one chunk's singles to a CSV text part-file via a file sink over `singles_chunk!`:
# format the row (positions cm→mm, energy MeV→keV) and stream it. Returns the row count.
function write_chunk_csv(path, geom, src, E0::Float64, cut_MeV::Float64, acol::Float64,
                         range::UnitRange{Int}, rng::AbstractRNG)
    open(path, "w") do io
        singles_chunk!(geom, src, E0, cut_MeV, acol, range, rng) do ev, g, s, pos0
            println(io, join((ev, g,
                round(s.x*10, digits=4), round(s.y*10, digits=4), round(s.z*10, digits=4),
                round(s.e*1000, digits=4), s.iz, s.iphi, s.nblocks, s.phscat ? 1 : 0,
                round(pos0[1]*10, digits=4), round(pos0[2]*10, digits=4), round(pos0[3]*10, digits=4)), ","))
        end
    end
end

# Write one chunk's singles to a quantized BINARY part-file: accumulate into a columnar
# SinglesBuffer (no HDF5 in the parallel region), then dump it. The serial glue packs the
# parts into HDF5. Returns the row count.
function write_chunk_bin(path, geom, src, E0::Float64, cut_MeV::Float64, acol::Float64,
                         range::UnitRange{Int}, rng::AbstractRNG)
    buf = SinglesBuffer()
    singles_chunk!(geom, src, E0, cut_MeV, acol, range, rng) do ev, g, s, pos0
        push_single!(buf, ev, g, s, pos0)
    end
    open(io -> write_part(io, buf), path, "w")
    length(buf)
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    datadir  = rp(cfg_get(cfg, "paths", "data", "data"))
    geomfile = rp(cfg_get(cfg, "geometry", "file", "geometry/geometry.json"))
    phmat    = String(cfg_get(cfg, "geometry", "phantom_material", ""))
    crystal  = String(cfg_get(cfg, "transport", "crystal_material", "CsI"))
    kind     = String(cfg_get(cfg, "source", "kind", "point"))
    acol     = Float64(cfg_get(cfg, "source", "acol_fwhm_deg", 0.5))
    energy   = Float64(cfg_get(cfg, "source", "energy_keV", 511.0))
    cutoff   = Float64(cfg_get(cfg, "transport", "cutoff_keV", 10.0))
    nevents  = a["nevents"] >= 0 ? a["nevents"] : Int(cfg_get(cfg, "source", "nevents", 100000))
    seed     = a["seed"]    >= 0 ? a["seed"]    : Int(cfg_get(cfg, "transport", "seed", 1234))
    fmt      = isempty(a["format"]) ? String(cfg_get(cfg, "output", "format", "csv")) : a["format"]

    kind in ("point", "phantom") || error("[source].kind must be 'point' or 'phantom'")
    fmt in ("csv", "hdf5") || error("[output].format '$fmt' not supported (csv | hdf5)")
    nevents >= 1 || error("nevents must be ≥ 1")
    ishdf5 = fmt == "hdf5"

    tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(prod_base(cfg)), tag)
    mkpath(outdir)
    cp(a["config"], joinpath(outdir, "config.toml"); force=true)
    out = joinpath(outdir, ishdf5 ? "singles.h5" : "singles.csv")

    mats = load_materials(datadir)
    geom = load_geometry(geomfile, mats)
    geom.scanner === nothing && error("geometry $geomfile has no scanner section")
    haskey(mats, crystal) || error("crystal material '$crystal' not found")

    # Crystal material swap, keeping the ring geometry/segmentation.
    sc0 = geom.scanner
    sc  = Scanner(PhysicalVolume(LogicalVolume(name(sc0.volume), solid(sc0.volume),
                                              mats[crystal]), sc0.volume.position),
                  sc0.n_phi, sc0.n_z)
    # Optional phantom-material override.
    phantom = geom.phantom
    if !isempty(phmat)
        haskey(mats, phmat) || error("phantom material '$phmat' not found")
        phantom = PhysicalVolume(LogicalVolume(name(phantom), solid(phantom), mats[phmat]),
                                 phantom.position)
    end
    geom = Geometry(geom.world, phantom, sc)

    E0 = energy / 1000.0
    cut_MeV = cutoff / 1000.0
    uniform = kind == "phantom"
    src = uniform ? UniformVolumeSource(geom.phantom) : PointSource(geom.phantom.position)

    nthr    = nthreads()
    ranges  = chunk_ranges(nevents, a["nchunks"] > 0 ? a["nchunks"] : 8 * nthr)
    nchunks = length(ranges)

    # Pre-allocate one RNG per chunk BEFORE the loop (chunk 1 = base seed); the hot loop only
    # draws from its chunk's RNG, never constructs one.
    rngs = [MersenneTwister(seed + (c - 1)) for c in 1:nchunks]
    partfile(c) = joinpath(outdir, ishdf5 ? "singles.part$(c).bin" : "singles.part$(c).csv")
    rows = zeros(Int, nchunks)

    srcdesc = uniform ? "uniform in $(name(geom.phantom)) ($(material(geom.phantom).name))" :
                        "point at the phantom centre ($(material(geom.phantom).name))"
    println("run '$tag': back-to-back $energy keV, $srcdesc, acol $(acol)° FWHM; " *
            "ring $(name(sc.volume)) ($crystal), nblocks=$(nblocks(sc)), cutoff $cutoff keV, " *
            "$nevents events, seed $seed; threads=$nthr, chunks=$nchunks, format=$fmt")

    write_chunk = ishdf5 ? write_chunk_bin : write_chunk_csv
    t0 = time()
    @threads :dynamic for c in 1:nchunks
        rows[c] = write_chunk(partfile(c), geom, src, E0, cut_MeV, acol, ranges[c], rngs[c])
    end
    t_sim = time() - t0

    # Glue the part-files in chunk order → one event-ordered output (bounded memory).
    nsingles = sum(rows)
    parts = [partfile(c) for c in 1:nchunks]
    if ishdf5
        meta = Dict{String,Any}("scenario_tag"=>tag, "seed"=>seed, "nevents"=>nevents,
            "nchunks"=>nchunks, "crystal"=>crystal, "phantom_material"=>material(geom.phantom).name,
            "energy_keV"=>energy, "acol_fwhm_deg"=>acol, "cutoff_keV"=>cutoff,
            "geometry_file"=>basename(geomfile), "n_phi"=>sc.n_phi, "n_z"=>sc.n_z)
        pack_singles_hdf5(out, parts, rows, meta)
    else
        open(out, "w") do io
            println(io, "event_number,gamma,x_mm,y_mm,z_mm,e_keV,iz,iphi,nblocks,phantom_scatter,x0_mm,y0_mm,z0_mm")
            buf = Vector{UInt8}(undef, 1 << 20)
            for p in parts
                open(p, "r") do pin
                    while !eof(pin)
                        n = readbytes!(pin, buf)
                        write(io, view(buf, 1:n))
                    end
                end
                rm(p)
            end
        end
    end

    ngamma = 2 * nevents
    dt = round(time() - t0, digits=2)
    println("wrote $nsingles singles ($(round(100 * nsingles / ngamma, digits=1))% of $ngamma photons) " *
            "-> $out")
    println("  transport $(round(t_sim, digits=2))s + glue $(round(time() - t0 - t_sim, digits=2))s = $(dt)s " *
            "($(round(nevents / max(t_sim, 1e-9) / 1000, digits=1))k events/s on $nthr threads)")
end

main()
