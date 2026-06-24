# Run configuration: the LOR-generation pipeline (simulate_phantom, build_coincidences)
# is driven by a TOML config — the parameter source of truth and the run's provenance.
# See docs/PTCryspMC_app.tex and runs/*.toml for the schema and directory layout.

import TOML

"Parse a run-config TOML file into a nested Dict."
read_config(path::AbstractString)::Dict = TOML.parsefile(path)

"""
    cfg_get(cfg, section, key, default)

Fetch `cfg[section][key]`, or `default` if the section or key is absent.
"""
cfg_get(cfg::AbstractDict, section::AbstractString, key::AbstractString, default) =
    get(get(cfg, section, Dict{String,Any}()), key, default)

"""
    run_tag(cfg, configpath) -> String

The run identifier: `[output].tag` if set, else the config filename without its extension.
Names the per-run output directory `output/<tag>/` and the stamped metadata.
"""
function run_tag(cfg::AbstractDict, configpath::AbstractString)::String
    t = cfg_get(cfg, "output", "tag", "")
    isempty(String(t)) ? splitext(basename(configpath))[1] : String(t)
end

"""
    prod_base(cfg) -> String

Base directory for the **production** chain outputs (singles + LORs), from `[output].prod_dir`
(default `"prod"`). Kept separate from the dev chain's `[output].dir` (default `"output"`, full
stacks + dev coincidences + plots) so the two pipelines never collide in one `output/<tag>/`.
"""
prod_base(cfg::AbstractDict)::String = String(cfg_get(cfg, "output", "prod_dir", "prod"))
