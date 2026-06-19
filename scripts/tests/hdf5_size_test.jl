#!/usr/bin/env julia
# Measurement: encode a singles CSV to the quantized-integer schema and write it as HDF5 at
# several compression settings, to see how much HDF5 saves over float/int CSV. Columnar typed
# datasets: event Int32, gamma/nblocks/phantom_scatter Int8, iz/iphi/positions/energy Int16
# (positions in 0.1 mm units, energy in 0.1 keV units; bounds-checked against Int16).
#
# Loads all columns into memory (fine for a ~10^7 run; this is an experiment, not the
# production writer). Config-driven; --in overrides the input path.
#   julia --project=. scripts/hdf5_size_test.jl --config runs/sphere_water_csi.toml

using PTCryspMC
using ArgParse
using HDF5

const SC_XYZ = 0.1; const SC_E = 0.1; const I16MAX = 32767

function parse_cli()
    s = ArgParseSettings(description="Measure HDF5 size of quantized-integer singles.")
    @add_arg_table! s begin
        "--config"; help = "run config TOML"; required = true
        "--in";     help = "override input singles path (default output/<tag>/singles.csv)"; default = ""
    end
    parse_args(s)
end

enc(v, sc) = (u = round(Int, v / sc); abs(u) <= I16MAX || error("overflow $v"); Int16(u))
mb(b) = "$(round(b/1_000_000, digits=1)) MB"

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))
    cfg = read_config(a["config"]); tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(prod_base(cfg)), tag)
    inp = isempty(a["in"]) ? joinpath(outdir, "singles.csv") : rp(a["in"])
    isfile(inp) || error("input '$inp' not found")

    ev = Int32[]; ga = Int8[]; x = Int16[]; y = Int16[]; z = Int16[]; e = Int16[]
    iz = Int16[]; ip = Int16[]; nb = Int8[]; ph = Int8[]; x0 = Int16[]; y0 = Int16[]; z0 = Int16[]
    open(inp) do io
        readline(io)
        for ln in eachline(io)
            isempty(ln) && continue
            f = split(ln, ',')
            push!(ev, parse(Int32, f[1])); push!(ga, parse(Int8, f[2]))
            push!(x, enc(parse(Float64, f[3]), SC_XYZ)); push!(y, enc(parse(Float64, f[4]), SC_XYZ)); push!(z, enc(parse(Float64, f[5]), SC_XYZ))
            push!(e, enc(parse(Float64, f[6]), SC_E))
            push!(iz, parse(Int16, f[7])); push!(ip, parse(Int16, f[8]))
            push!(nb, parse(Int8, f[9])); push!(ph, parse(Int8, f[10]))
            push!(x0, enc(parse(Float64, f[11]), SC_XYZ)); push!(y0, enc(parse(Float64, f[12]), SC_XYZ)); push!(z0, enc(parse(Float64, f[13]), SC_XYZ))
        end
    end
    n = length(ev); println("rows = $n   (encoded, no Int16 overflow)")

    cols = [("event_number", ev), ("gamma", ga), ("x_mm", x), ("y_mm", y), ("z_mm", z),
            ("e_keV", e), ("iz", iz), ("iphi", ip), ("nblocks", nb), ("phantom_scatter", ph),
            ("x0_mm", x0), ("y0_mm", y0), ("z0_mm", z0)]
    C = 1 << 16
    variants = (("uncompressed", (; chunk=(C,))),
                ("deflate=4",    (; chunk=(C,), deflate=4)),
                ("deflate=9",    (; chunk=(C,), deflate=9)),
                ("shuffle+deflate=4", (; chunk=(C,), shuffle=true, deflate=4)),
                ("shuffle+deflate=9", (; chunk=(C,), shuffle=true, deflate=9)))
    for (name, opts) in variants
        p = joinpath(outdir, "singles_sizetest.h5")
        t = @elapsed h5open(p, "w") do h
            for (nm, v) in cols
                d = create_dataset(h, nm, eltype(v), (length(v),); opts...)
                write(d, v)
            end
        end
        println("  HDF5 $name: $(mb(filesize(p)))   (write $(round(t, digits=1))s)")
        rm(p)
    end
    println("reference: float CSV $(mb(filesize(inp)))   (int CSV gzip-6 ≈ 256 MB)")
end

main()
