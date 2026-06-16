#!/usr/bin/env julia
# Phantom simulation driver — TOML-config driven. A back-to-back 511 keV source fills the
# phantom (uniform by volume, or a point at its centre) and emits pairs with acollinearity;
# each photon is NAVIGATED across the full geometry — phantom → air gap → crystal ring — by
# `navigate_photon`. We record the full interaction stack of both photons, tagged with the
# volume (:phantom / :air / :scanner), the crystal (iz, iphi) of each scanner deposit, a
# per-photon phantom-scatter flag, and the true emission point (x0, y0, z0).
#
# The phantom is a solid (cylinder or sphere) + material from the geometry JSON, with the
# material optionally overridden (Vacuum/Air = non-attenuated reference, Water = realistic).
# Source = attenuator: the same volume both emits and attenuates.
#
# All parameters come from a run config TOML (see dev/phantom_track_plan.md, runs/*.toml);
# the run `tag` (= config filename) names the per-run output dir output/<tag>/, where the
# stack and a copy of the config are written. Its air-only counterpart is shoot_into_ring.jl.
#
# Run from the repo root:
#   julia --project=. scripts/simulate_phantom.jl --config runs/sphere_water_csi.toml
#   julia --project=. scripts/simulate_phantom.jl --config runs/sphere_vacuum_csi.toml --nevents 5000

using PTCryspMC
using ArgParse
using Random

function parse_cli()
    s = ArgParseSettings(description="Phantom simulation driver (TOML-config driven).")
    @add_arg_table! s begin
        "--config";  help = "run config TOML"; required = true
        "--nevents"; help = "override [source].nevents"; arg_type = Int; default = -1
        "--seed";    help = "override [transport].seed"; arg_type = Int; default = -1
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
    # TRUTH full-energy: ≥505 keV in one crystal (implies ~unscattered in the phantom).
    full_E   = Escan >= 0.505 && length(blocks) == 1
    (reached, full_E, phscat, Escan, length(blocks))
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    datadir = rp(cfg_get(cfg, "paths", "data", "data"))
    geomfile = rp(cfg_get(cfg, "geometry", "file", "geometry/geometry.json"))
    phmat    = String(cfg_get(cfg, "geometry", "phantom_material", ""))
    crystal  = String(cfg_get(cfg, "transport", "crystal_material", "CsI"))
    kind     = String(cfg_get(cfg, "source", "kind", "point"))
    acol     = Float64(cfg_get(cfg, "source", "acol_fwhm_deg", 0.5))
    energy   = Float64(cfg_get(cfg, "source", "energy_keV", 511.0))
    cutoff   = Float64(cfg_get(cfg, "transport", "cutoff_keV", 10.0))
    nevents  = a["nevents"] >= 0 ? a["nevents"] : Int(cfg_get(cfg, "source", "nevents", 100000))
    seed     = a["seed"]    >= 0 ? a["seed"]    : Int(cfg_get(cfg, "transport", "seed", 1234))
    fmt      = String(cfg_get(cfg, "output", "format", "csv"))
    singles  = Bool(cfg_get(cfg, "output", "singles", false))

    kind in ("point", "phantom") || error("[source].kind must be 'point' or 'phantom'")
    fmt == "csv" || error("[output].format '$fmt' not supported yet (csv only; hdf5 deferred)")
    singles && error("[output].singles not implemented yet (full stack only; deferred)")

    tag = run_tag(cfg, a["config"])
    outdir = joinpath(rp(cfg_get(cfg, "output", "dir", "output")), tag)
    mkpath(outdir)
    cp(a["config"], joinpath(outdir, "config.toml"); force=true)
    out = joinpath(outdir, "stack.csv")

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
    rng = MersenneTwister(seed)
    src = uniform ? UniformVolumeSource(geom.phantom) : PointSource(geom.phantom.position)

    srcdesc = uniform ? "uniform in $(name(geom.phantom)) ($(material(geom.phantom).name))" :
                        "point at the phantom centre ($(material(geom.phantom).name))"
    println("run '$tag': back-to-back $energy keV, $srcdesc, acol $(acol)° FWHM; " *
            "ring $(name(sc.volume)) ($crystal), nblocks=$(nblocks(sc)), cutoff $cutoff keV, " *
            "$nevents events, seed $seed")

    n_reach = 0; n_fullE = 0; n_phscat = 0
    n_both_fullE = 0; n_both_unscat = 0
    nrows = 0
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
                fe = g == 1 ? (full_E, fe[2])  : (fe[1], full_E)
                us = g == 1 ? (!phscat, us[2]) : (us[1], !phscat)
                for (k, st) in enumerate(steps)
                    r = st.hit
                    println(io, join((ev, g, k,
                        round(r.x*10, digits=4), round(r.y*10, digits=4), round(r.z*10, digits=4),
                        round(r.e_in*1000, digits=4), round(r.e_dep*1000, digits=4), r.process,
                        st.volume, st.iz, st.iphi, phscat ? 1 : 0, x0, y0, z0), ","))
                    nrows += 1
                end
            end
            n_both_unscat += (us[1] && us[2])
            n_both_fullE  += (fe[1] && fe[2])
        end
    end

    ngamma = 2 * nevents
    println("wrote $nevents events ($nrows interaction rows) -> $out")
    println("per-photon:  reached ring $(round(100*n_reach/ngamma, digits=1))%   " *
            "unscattered in phantom $(round(100*(1 - n_phscat/ngamma), digits=1))%   " *
            "full-energy (≥505 keV, 1 crystal) $(round(100*n_fullE/ngamma, digits=1))%")
    println("per-event:   both photons unscattered $(round(100*n_both_unscat/nevents, digits=1))%   " *
            "both full-energy (truth) $(round(100*n_both_fullE/nevents, digits=1))%")
    println("NOTE: 'full-energy' is a TRUTH ≥505 keV cut — the UNSCATTERED, fully-contained")
    println("      subset, NOT the coincidence efficiency (that comes from build_coincidences).")
end

main()
