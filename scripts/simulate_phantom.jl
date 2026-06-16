#!/usr/bin/env julia
# Phantom simulation driver. A back-to-back 511 keV source fills the phantom (uniform by
# volume, or a point at its centre) and emits pairs with ~0.5 deg FWHM acollinearity; each
# photon is NAVIGATED across the full geometry — phantom → air gap → crystal ring — by
# `navigate_photon`. We record the full interaction stack of both photons, tagged with the
# volume (:phantom / :air / :scanner), the crystal (iz, iphi) of each scanner deposit, a
# per-photon phantom-scatter flag, and the true emission point (x0, y0, z0).
#
# The phantom is a solid (cylinder or sphere) + material: --phantom-material picks it, so
# Vacuum/Air is the non-attenuated reference (no scatter) and Water is realistic. Source =
# attenuator: the same volume both emits and attenuates (see dev/phantom_track_plan.md).
# Its air-only counterpart is shoot_into_ring.jl (point source, no phantom).
#
# IMPORTANT — what "full-energy" means here. These are TRUTH energies (no detector
# resolution yet). A photon is called full-energy if it deposits ≥505 keV in a single
# crystal. At truth level a fully-absorbed photon deposits exactly 511 keV, so this just
# picks out the full-energy peak — BUT a ≥505 keV cut also rejects every photon that
# Compton-scattered in the phantom (it arrives with <511 keV). So the "both full-energy"
# fraction below is the UNSCATTERED, fully-contained subset — NOT the coincidence
# efficiency. The real selection (Step 4) smears by FWHM(E) and applies a ±2·FWHM energy
# window, which also accepts scattered photons as *scatter* coincidences.
#
# Run from the repo root:
#   julia --project=. scripts/simulate_phantom.jl --source phantom --material CsI
#   julia --project=. scripts/simulate_phantom.jl --geometry geometry/geometry_sphere.json \
#         --source phantom --phantom-material Vacuum --material BGO   # non-attenuated reference

using PTCryspMC
using ArgParse
using Random

function parse_cli()
    s = ArgParseSettings(description="Back-to-back 511 keV pairs from the phantom, navigated through water → air → ring.")
    @add_arg_table! s begin
        "--data";     help = "data dir";        default = joinpath(@__DIR__, "..", "data")
        "--geometry"; help = "geometry JSON";   default = joinpath(@__DIR__, "..", "geometry", "geometry.json")
        "--out";      help = "output CSV (default output/phantom_<material>_stack.csv)"; default = ""
        "--material"; help = "crystal material"; default = "CsI"
        "--phantom-material"; help = "override the phantom material (default = geometry JSON value; Vacuum = non-attenuated reference)"; default = ""
        "--source";   help = "source: 'point' (phantom centre) or 'phantom' (uniform in the phantom)"; default = "point"
        "--acol-fwhm"; help = "acollinearity FWHM [deg] (0 = exactly back-to-back)"; arg_type = Float64; default = 0.5
        "--nevents";  help = "n annihilations"; arg_type = Int;     default = 100000
        "--energy";   help = "energy [keV]";    arg_type = Float64; default = 511.0
        "--cutoff";   help = "low-energy cutoff [keV] (default 10 = XCOM min)"; arg_type = Float64; default = 10.0
        "--seed";     help = "RNG seed";        arg_type = Int;     default = 1234
    end
    parse_args(s)
end

# Per-photon summary derived from its navigated stack.
function photon_stats(steps)
    Escan = 0.0; phscat = false
    blocks = Set{Tuple{Int,Int}}()
    for st in steps
        if st.volume == :phantom && st.hit.process != :escape
            phscat = true
        elseif st.volume == :scanner && st.hit.process != :escape
            Escan += st.hit.e_dep
            push!(blocks, (st.iz, st.iphi))
        end
    end
    reached  = !isempty(blocks)
    # TRUTH full-energy: ≥505 keV in one crystal. At truth level a fully-absorbed photon
    # deposits exactly 511 keV, so this is the full-energy peak; it also implies the photon
    # was ~unscattered in the phantom (a phantom Compton would leave it below 505).
    full_E   = Escan >= 0.505 && length(blocks) == 1
    (reached, full_E, phscat, Escan, length(blocks))
end

function main()
    a = parse_cli()
    mats = load_materials(a["data"])
    geom = load_geometry(a["geometry"], mats)
    geom.scanner === nothing && error("geometry $(a["geometry"]) has no scanner section")
    haskey(mats, a["material"]) || error("material '$(a["material"])' not found")
    a["source"] in ("point", "phantom") || error("--source must be 'point' or 'phantom'")

    # Swap in the chosen crystal material, keeping the ring geometry/segmentation.
    sc0 = geom.scanner
    sc  = Scanner(PhysicalVolume(LogicalVolume(name(sc0.volume), solid(sc0.volume),
                                               mats[a["material"]]), sc0.volume.position),
                  sc0.n_phi, sc0.n_z)

    # Optional phantom-material override (Vacuum = non-attenuated reference, Water = realistic).
    phantom = geom.phantom
    if !isempty(a["phantom-material"])
        pm = a["phantom-material"]
        haskey(mats, pm) || error("phantom material '$pm' not found")
        phantom = PhysicalVolume(LogicalVolume(name(phantom), solid(phantom), mats[pm]),
                                 phantom.position)
    end
    geom = Geometry(geom.world, phantom, sc)

    out = isempty(a["out"]) ?
        joinpath(@__DIR__, "..", "output", "phantom_$(lowercase(a["material"]))_stack.csv") : a["out"]
    E0      = a["energy"] / 1000.0
    cut_MeV = a["cutoff"] / 1000.0
    uniform = a["source"] == "phantom"
    acol    = a["acol-fwhm"]
    rng     = MersenneTwister(a["seed"])

    # Source = attenuator: fill the phantom uniformly (or a point at its centre).
    src = uniform ? UniformVolumeSource(geom.phantom) : PointSource(geom.phantom.position)

    srcdesc = uniform ? "uniform in $(name(geom.phantom)) ($(material(geom.phantom).name))" :
                        "point at the phantom centre ($(material(geom.phantom).name))"
    println("source: back-to-back $(a["energy"]) keV, $srcdesc, acol $(acol)° FWHM; " *
            "ring $(name(sc.volume)) ($(a["material"])), nblocks=$(nblocks(sc)), " *
            "cutoff $(a["cutoff"]) keV")

    n_reach = 0; n_fullE = 0; n_phscat = 0
    n_both_fullE = 0; n_both_unscat = 0
    nevents = a["nevents"]
    nrows = 0
    mkpath(dirname(out))
    open(out, "w") do io
        println(io, "event_number,gamma,step,x_mm,y_mm,z_mm,e_in_keV,e_dep_keV,process,volume,iz,iphi,phantom_scatter,x0_mm,y0_mm,z0_mm")
        for ev in 1:nevents
            pos0, d1, d2 = emit_pair(src, rng; acol_fwhm_deg=acol)
            x0 = round(pos0[1]*10, digits=4); y0 = round(pos0[2]*10, digits=4); z0 = round(pos0[3]*10, digits=4)
            fe = (false, false); us = (false, false)
            for (g, dir) in ((1, d1), (2, d2))
                steps = navigate_photon(geom, E0, pos0, dir, rng; egamma_cut=cut_MeV)
                reached, full_E, phscat, _, _ = photon_stats(steps)
                n_reach  += reached
                n_fullE  += full_E
                n_phscat += phscat
                fe = g == 1 ? (full_E, fe[2])      : (fe[1], full_E)
                us = g == 1 ? (!phscat, us[2])     : (us[1], !phscat)
                for (k, st) in enumerate(steps)
                    r = st.hit
                    println(io, join((ev, g, k,
                        round(r.x*10, digits=4), round(r.y*10, digits=4), round(r.z*10, digits=4),
                        round(r.e_in*1000, digits=4), round(r.e_dep*1000, digits=4), r.process,
                        st.volume, st.iz, st.iphi, phscat ? 1 : 0, x0, y0, z0), ","))
                    nrows += 1
                end
            end
            n_both_unscat += (us[1] && us[2])      # both photons crossed the phantom unscattered
            n_both_fullE  += (fe[1] && fe[2])      # both fully contained at 511 keV (truth)
        end
    end

    ngamma = 2 * nevents
    println("wrote $nevents events ($nrows interaction rows) -> $out")
    println("per-photon:  reached ring $(round(100*n_reach/ngamma, digits=1))%   " *
            "unscattered in phantom $(round(100*(1 - n_phscat/ngamma), digits=1))%   " *
            "full-energy (≥505 keV, 1 crystal) $(round(100*n_fullE/ngamma, digits=1))%")
    println("per-event:   both photons unscattered $(round(100*n_both_unscat/nevents, digits=1))%   " *
            "both full-energy (truth) $(round(100*n_both_fullE/nevents, digits=1))%  " *
            "($n_both_fullE / $nevents)")
    println("NOTE: 'full-energy' is a TRUTH ≥505 keV cut — it is the UNSCATTERED, fully-")
    println("      contained subset, NOT the coincidence efficiency. The real ±2·FWHM energy")
    println("      window (Step 4) also keeps scattered photons as *scatter* coincidences.")
end

main()
