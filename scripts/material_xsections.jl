#!/usr/bin/env julia
# Dump the macroscopic photon cross sections (Compton, photoelectric, pair) of a
# material over a log-spaced energy grid, for a companion Python plot script.
# The sampling is decoupled from the XCOM grid: sigma_macro is interpolated at
# `npoints` energies between `emin` and `emax`.
#
# Run from the repo root:
#   julia --project=. scripts/material_xsections.jl --material Water
#   julia --project=. scripts/material_xsections.jl --material CsI
#
# With no --out, the file is output/<material lowercased>_xsections.csv.

using PTCryspMC
using ArgParse
using Printf

function parse_cli()
    s = ArgParseSettings(description="Tabulate macroscopic photon cross sections [cm^-1] vs energy for a material.")
    @add_arg_table! s begin
        "--data";     help = "data dir";          default = joinpath(@__DIR__, "..", "data")
        "--material"; help = "material name";     default = "Water"
        "--out";      help = "output CSV (default output/<material>_xsections.csv)"; default = ""
        "--emin";     help = "min energy [keV]";  arg_type = Float64; default = 10.0
        "--emax";     help = "max energy [keV]";  arg_type = Float64; default = 10000.0
        "--npoints";  help = "n grid points";     arg_type = Int;     default = 20
    end
    parse_args(s)
end

function main()
    a = parse_cli()
    mat = load_material(a["data"], a["material"])
    out = isempty(a["out"]) ?
        joinpath(@__DIR__, "..", "output", lowercase(a["material"]) * "_xsections.csv") :
        a["out"]

    emin_MeV, emax_MeV = a["emin"] / 1000.0, a["emax"] / 1000.0
    # The table is only valid inside its tabulated range; refuse to extrapolate.
    lo, hi = first(mat.E), last(mat.E)
    (emin_MeV >= lo && emax_MeV <= hi) ||
        error("requested range [$(a["emin"]), $(a["emax"])] keV is outside the " *
              "$(mat.name) XCOM table [$(lo*1000), $(hi*1000)] keV")

    # Log-spaced grid (cross sections span decades; a log axis is the natural plot).
    n = a["npoints"]
    logE = range(log10(emin_MeV), log10(emax_MeV); length = n)

    mkpath(dirname(out))
    open(out, "w") do io
        println(io, "energy_keV,compton,phot,pair")   # all sigma in cm^-1 (macroscopic)
        for le in logE
            E = 10.0^le
            ΣC, ΣPh, ΣP = sigma_macro(mat, E)
            @printf(io, "%.6g,%.6g,%.6g,%.6g\n", E * 1000, ΣC, ΣPh, ΣP)
        end
    end
    println("wrote $n points for $(mat.name), $(a["emin"])–$(a["emax"]) keV -> $out")
end

main()
