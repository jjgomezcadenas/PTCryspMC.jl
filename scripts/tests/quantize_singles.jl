#!/usr/bin/env julia
# Quantize a singles CSV to INTEGER columns and report the size reduction (the float→int
# experiment, made reproducible). Positions are stored in 0.1 mm units and energy in 0.1 keV
# units as integers — far finer than the ~2 mm / ~50 keV detector resolution, so the
# quantization is lossless in practice (and downstream smearing dwarfs the grid). Each value
# is bounds-checked against the Int16 range (±3276.7 mm / 0–3276.7 keV) so an oversized
# geometry fails loud instead of silently wrapping.
#
# Streams in / out (O(1) memory) and reports raw + gzip sizes for the float input and the int
# output, so the on-disk saving is measured, not estimated.
#
# TOML-config driven (reads output/<tag>/singles.csv), mirroring check_singles.jl:
#   julia --project=. scripts/quantize_singles.jl --config runs/sphere_water_csi.toml
#   julia --project=. scripts/quantize_singles.jl --config runs/sphere_water_csi.toml --in path.csv --out q.csv

using PTCryspMC
using ArgParse

const XYZ_SCALE_MM = 0.1     # mm per integer unit  → int = round(mm / 0.1)
const E_SCALE_KEV  = 0.1     # keV per integer unit
const I16MAX = 32767

# Columns quantized (mm / keV → scaled integer) vs passed through unchanged (already integer).
const QUANT_COLS = ("x_mm", "y_mm", "z_mm", "e_keV", "x0_mm", "y0_mm", "z0_mm")
const SCALE_OF   = Dict("x_mm"=>XYZ_SCALE_MM, "y_mm"=>XYZ_SCALE_MM, "z_mm"=>XYZ_SCALE_MM,
                        "e_keV"=>E_SCALE_KEV, "x0_mm"=>XYZ_SCALE_MM, "y0_mm"=>XYZ_SCALE_MM,
                        "z0_mm"=>XYZ_SCALE_MM)

function parse_cli()
    s = ArgParseSettings(description="Quantize a singles CSV to integers and report the size saving.")
    @add_arg_table! s begin
        "--config"; help = "run config TOML"; required = true
        "--in";     help = "override the input singles path (default output/<tag>/singles.csv)"; default = ""
        "--out";    help = "override the output integer-CSV path (default <in dir>/singles_int.csv)"; default = ""
    end
    parse_args(s)
end

# Gzip a file via the system gzip and return the compressed byte count (-1 if gzip absent).
function gzip_size(path, level)
    try
        parse(Int, strip(read(pipeline(`gzip -$level -c $path`, `wc -c`), String)))
    catch
        -1
    end
end

mb(b) = b < 0 ? "n/a" : "$(round(b/1_000_000, digits=1)) MB"

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(cfg_get(cfg, "output", "dir", "output")), tag)
    inp = isempty(a["in"])  ? joinpath(outdir, "singles.csv") : rp(a["in"])
    isfile(inp) || error("input '$inp' not found (run simulate_source_mt.jl first)")
    out = isempty(a["out"]) ? joinpath(dirname(inp), "singles_int.csv") : rp(a["out"])

    rows = 0
    maxabs_mm = 0.0; maxe_keV = 0.0
    open(inp, "r") do fin
        header = split(strip(readline(fin)), ',')
        col = Dict(String(h) => i for (i, h) in enumerate(header))
        for c in QUANT_COLS
            haskey(col, c) || error("singles stack is missing column '$c'")
        end
        qidx = Set(col[c] for c in QUANT_COLS)                  # which 1-based fields to scale

        open(out, "w") do fout
            println(fout, join(header, ","))                     # same column names
            for line in eachline(fin)
                isempty(line) && continue
                f = split(line, ',')
                rows += 1
                for i in eachindex(f)
                    i == 1 || print(fout, ',')
                    if i in qidx
                        v = parse(Float64, f[i])
                        sc = SCALE_OF[String(header[i])]
                        u = round(Int, v / sc)
                        abs(u) <= I16MAX || error("value $v in column '$(header[i])' exceeds " *
                            "Int16 @ $sc (±$(I16MAX*sc)) at row $rows — widen the scale")
                        print(fout, u)
                        a_mm = abs(v); a_mm > maxabs_mm && header[i] != "e_keV" && (maxabs_mm = a_mm)
                        header[i] == "e_keV" && v > maxe_keV && (maxe_keV = v)
                    else
                        print(fout, f[i])                        # passthrough (already integer)
                    end
                end
                println(fout)
            end
        end
    end

    fsz = filesize(inp); isz = filesize(out)
    println("quantized $rows rows: $inp")
    println("  -> $out")
    println("  range seen: |pos| ≤ $(round(maxabs_mm, digits=1)) mm, E ≤ $(round(maxe_keV, digits=1)) keV " *
            "(Int16 @0.1 covers ±$(I16MAX*XYZ_SCALE_MM) mm / $(I16MAX*E_SCALE_KEV) keV) — no overflow")
    println("  raw:      float $(mb(fsz))   int $(mb(isz))   ($(round(100*(1-isz/fsz), digits=0))% smaller)")
    gf1 = gzip_size(inp, 1); gi1 = gzip_size(out, 1)
    gf6 = gzip_size(inp, 6); gi6 = gzip_size(out, 6)
    println("  gzip -1:  float $(mb(gf1))   int $(mb(gi1))")
    println("  gzip -6:  float $(mb(gf6))   int $(mb(gi6))")
end

main()
