#!/usr/bin/env julia
# Publish a completed API production run into the PtCryspProds/ products tree (the downstream
# handoff). Derives the layout leaf from the config —
#     <scenario>/<scanner>/<crystal>/<budget>_<dose>/lors_shard<NNN>.h5
# — copies the run's lors_det.h5 there as a shard, and writes the shared files once: config.toml at
# the leaf, scanner_geometry.json at the scanner level, phantom/ (the μ-map inputs) at the scenario
# level, and README.md (the contract) at the root. Idempotent; --force overwrites an existing shard.
# See dev/PRODUCTS.md for the layout contract and dev/data_generation_strategy.md for the meaning.
#
#   julia --project=. scripts/run/publish_prod.jl --config runs/uniform_headep_bgo_api.toml \
#         --prods-root ~/Projects/PtCryspProds

using PTCryspMC
using ArgParse
import JSON

function parse_cli()
    s = ArgParseSettings(description="Publish an API prod run into the PtCryspProds tree.")
    @add_arg_table! s begin
        "--config";      help = "the API run config TOML (mode=api)"; required = true
        "--prods-root";  help = "the PtCryspProds root directory"; required = true
        "--force";       help = "overwrite an existing shard file"; action = :store_true
        "--truth-only";  help = "write only the scenario truth/ bundle (backfill; no shard/lors needed)"; action = :store_true
        "--force-truth"; help = "overwrite an existing truth/ bundle"; action = :store_true
    end
    parse_args(s)
end

# dose_Gy → a filesystem-safe token: 1.0→"1Gy", 0.5→"0p5Gy", 2.5→"2p5Gy".
dose_tag(d) = (s = d == round(d) ? string(Int(round(d))) : replace(string(d), "." => "p"); s * "Gy")
sanitize(s) = lowercase(replace(String(s), r"[^0-9A-Za-z]+" => "_"))

# Copy `src`→`dst` only if `dst` is absent (shared files written once); report what happened.
function copy_once(src, dst; label)
    if isfile(dst)
        println("    = $label (exists, kept)")
    else
        cp(src, dst; force=false); println("    + $label")
    end
end

function main()
    a = parse_cli()
    REPO = normpath(joinpath(@__DIR__, "..", ".."))
    rp(p) = (q = String(p); isabspath(q) ? q : joinpath(REPO, q))

    cfg = read_config(a["config"])
    String(cfg_get(cfg, "source", "mode", "")) == "api" || error("[source].mode must be 'api'")
    tag    = run_tag(cfg, a["config"])
    outdir = joinpath(rp(prod_base(cfg)), tag)

    # Layout components from the config (lors_det.h5 also carries these as provenance).
    scndir    = rp(String(cfg_get(cfg, "source", "scenario_dir", "")))
    scenario  = basename(rstrip(scndir, '/'))
    crystal   = sanitize(cfg_get(cfg, "transport", "crystal_material", "CsI"))
    budget    = String(cfg_get(cfg, "source", "budget", "fast"))
    dose_Gy   = Float64(cfg_get(cfg, "source", "dose_Gy", 1.0))
    dose      = dose_tag(dose_Gy)

    root      = rp(a["prods-root"])
    scen_dir  = joinpath(root, scenario)
    truth     = joinpath(scen_dir, "truth")

    # truth/ bundle: detector-independent, <scenario>/ level, shared like phantom/. Written once
    # per scenario (kept if present, unless --force-truth). Streams emitters.csv for activity(z).
    function publish_truth()
        mats = load_materials(rp(cfg_get(cfg, "paths", "data", "data")))
        println("  truth/ (dose-R80 + activity-R50, budget=$budget, dose=$(dose)):")
        before = isdir(truth) ? Set(readdir(truth)) : Set{String}()
        write_truth_bundle(scndir, truth, mats; budget=budget, dose_Gy=dose_Gy, force=a["force-truth"])
        for f in sort(readdir(truth))
            println(f in before && !a["force-truth"] ? "    = truth/$f (exists, kept)" : "    + truth/$f")
        end
    end

    if a["truth-only"]
        println("publishing truth-only for '$scenario' → $root")
        mkpath(scen_dir); publish_truth(); println("  done.")
        return
    end

    lors   = joinpath(outdir, "lors_det.h5")
    isfile(lors) || error("missing $lors — run (and don't prune) the chain first")
    realiz   = Int(singles_hdf5_attr(lors, "realization", cfg_get(cfg, "source", "realization", 0)))  # what was actually run
    geomfile = rp(cfg_get(cfg, "geometry", "file", "geometry/geometry_head.json"))
    scanner  = sanitize(get(get(JSON.parsefile(geomfile), "scanner", Dict()), "name", "scanner"))

    scan_dir  = joinpath(scen_dir, scanner)
    leaf      = joinpath(scan_dir, crystal, "$(budget)_$(dose)")
    phantom   = joinpath(scen_dir, "phantom")
    shard     = joinpath(leaf, "lors_shard$(lpad(realiz, 3, '0')).h5")

    println("publishing '$tag' → $root")
    println("  leaf: $scenario/$scanner/$crystal/$(budget)_$(dose)/  (shard $(lpad(realiz,3,'0')))")
    mkpath(leaf); mkpath(phantom)

    # The shard (the deliverable): copy lors_det.h5 → lors_shard<NNN>.h5.
    if isfile(shard) && !a["force"]
        println("    ! shard exists: $(basename(shard)) — use --force to overwrite (SKIPPED)")
    else
        cp(lors, shard; force=true); println("    + $(basename(shard))  ($(round(filesize(shard)/1e6,digits=1)) MB)")
    end

    # config.toml at the leaf (the recipe; shards differ only in realization).
    copy_once(joinpath(outdir, "config.toml"), joinpath(leaf, "config.toml"); label="config.toml")

    # scanner_geometry.json at the scanner level (the system model; shared across crystals).
    copy_once(geomfile, joinpath(scan_dir, "scanner_geometry.json"); label="../scanner_geometry.json")

    # phantom/ at the scenario level (the μ-map inputs; shared across scanners/crystals).
    for f in readdir(scndir)
        (f == "phantom_regions.csv" || startswith(f, "phantom_material_")) &&
            copy_once(joinpath(scndir, f), joinpath(phantom, f); label="../../phantom/$f")
    end

    # truth/ at the scenario level (dose + activity truth; shared across scanners/crystals).
    publish_truth()

    # README.md at the root (the contract), from dev/PRODUCTS.md.
    doc = rp("dev/PRODUCTS.md")
    isfile(doc) && copy_once(doc, joinpath(root, "README.md"); label="README.md (contract)")

    println("  done.")
end

main()
