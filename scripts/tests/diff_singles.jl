#!/usr/bin/env julia
# Compare two singles stacks for DATA equality — the reproducibility / regression check (e.g.
# -t 1 vs -t 18 at fixed nchunks, or CSV vs HDF5). HDF5 files can't be byte-diffed (internal
# metadata differs), so this compares the decoded columns: integer columns exactly, the
# physical columns (positions mm, energy keV) within `--tol`. Format is auto-detected per file
# from the extension (.csv / .h5), so any combination works (same-format runs match exactly;
# CSV-vs-HDF5 needs --tol ≈ 0.05 for the 0.1 unit quantization). Exits nonzero on any mismatch.
#
#   julia --project=. scripts/diff_singles.jl a/singles.h5 b/singles.h5
#   julia --project=. scripts/diff_singles.jl run/singles.csv run/singles.h5 --tol 0.05
#
# Loads both files into memory (decoded) — meant for validation-scale runs, not 10^8.

using PTCryspMC
using ArgParse

function parse_cli()
    s = ArgParseSettings(description="Compare two singles stacks (CSV/HDF5) for data equality.")
    @add_arg_table! s begin
        "a";     help = "first singles file (.csv or .h5)";  required = true
        "b";     help = "second singles file (.csv or .h5)"; required = true
        "--tol"; help = "tolerance [mm/keV] for the physical columns"; arg_type = Float64; default = 0.0
    end
    parse_args(s)
end

# Load a singles file into decoded columns: integer fields as Int, positions [mm] and energy
# [keV] as Float64 (HDF5 ints decoded via the stored scale).
function load_singles(path)
    ev = Int[]; ga = Int[]; x = Float64[]; y = Float64[]; z = Float64[]; e = Float64[]
    iz = Int[]; ip = Int[]; nb = Int[]; ph = Int[]; x0 = Float64[]; y0 = Float64[]; z0 = Float64[]
    if endswith(path, ".h5")
        foreach_singles_hdf5(path) do b
            for i in 1:length(b)
                push!(ev, b.event[i]); push!(ga, b.gamma[i])
                push!(x, decode_xyz(b.x[i])); push!(y, decode_xyz(b.y[i])); push!(z, decode_xyz(b.z[i]))
                push!(e, decode_e(b.e[i]))
                push!(iz, b.iz[i]); push!(ip, b.iphi[i]); push!(nb, b.nblocks[i]); push!(ph, b.phantom_scatter[i])
                push!(x0, decode_xyz(b.x0[i])); push!(y0, decode_xyz(b.y0[i])); push!(z0, decode_xyz(b.z0[i]))
            end
        end
    else
        open(path) do io
            header = split(strip(readline(io)), ',')
            c = Dict(String(h) => i for (i, h) in enumerate(header))
            for line in eachline(io)
                isempty(line) && continue
                f = split(line, ',')
                push!(ev, parse(Int, f[c["event_number"]])); push!(ga, parse(Int, f[c["gamma"]]))
                push!(x, parse(Float64, f[c["x_mm"]])); push!(y, parse(Float64, f[c["y_mm"]])); push!(z, parse(Float64, f[c["z_mm"]]))
                push!(e, parse(Float64, f[c["e_keV"]]))
                push!(iz, parse(Int, f[c["iz"]])); push!(ip, parse(Int, f[c["iphi"]]))
                push!(nb, parse(Int, f[c["nblocks"]])); push!(ph, parse(Int, f[c["phantom_scatter"]]))
                push!(x0, parse(Float64, f[c["x0_mm"]])); push!(y0, parse(Float64, f[c["y0_mm"]])); push!(z0, parse(Float64, f[c["z0_mm"]]))
            end
        end
    end
    (event=ev, gamma=ga, x=x, y=y, z=z, e=e, iz=iz, iphi=ip, nblocks=nb, phantom_scatter=ph, x0=x0, y0=y0, z0=z0)
end

function main()
    a = parse_cli()
    rp(p) = (q = String(p); isabspath(q) ? q : abspath(q))
    pa = rp(a["a"]); pb = rp(a["b"]); tol = a["tol"]
    isfile(pa) || error("file '$pa' not found"); isfile(pb) || error("file '$pb' not found")

    A = load_singles(pa); B = load_singles(pb)
    println("a: $pa  ($(length(A.event)) rows)")
    println("b: $pb  ($(length(B.event)) rows)")
    if length(A.event) != length(B.event)
        println("DIFFER ✗  (row counts $(length(A.event)) ≠ $(length(B.event)))"); exit(1)
    end

    intcols = (:event, :gamma, :iz, :iphi, :nblocks, :phantom_scatter)
    fltcols = (:x, :y, :z, :e, :x0, :y0, :z0)
    nbad = 0
    for col in intcols
        va = getfield(A, col); vb = getfield(B, col)
        m = findfirst(i -> va[i] != vb[i], eachindex(va))
        m === nothing || (nbad += 1; println("  col $col: first mismatch at row $m ($(va[m]) ≠ $(vb[m]))"))
    end
    for col in fltcols
        va = getfield(A, col); vb = getfield(B, col)
        m = findfirst(i -> abs(va[i] - vb[i]) > tol, eachindex(va))
        m === nothing || (nbad += 1; println("  col $col: first |Δ|>$tol at row $m ($(va[m]) vs $(vb[m]))"))
    end

    if nbad == 0
        println("IDENTICAL ✓  ($(length(A.event)) rows, tol=$tol)")
    else
        println("DIFFER ✗  ($nbad column(s) mismatched)"); exit(1)
    end
end

main()
