#!/usr/bin/env julia
# Generate docs/SCHEMA.md from the code, so the schema doc and the code cannot drift. The column
# names and types are introspected from the live schema (`singles_columns` / `coinc_columns` and
# the struct field eltypes) — NOT parsed from source text; the per-column units/meaning come from
# the co-located `singles_doc` / `coinc_doc` maps. A test (test/runtests.jl) regenerates this in
# memory and fails if docs/SCHEMA.md is stale or a column is undocumented.
#
#   julia --project=. scripts/gen_schema.jl        # (re)write docs/SCHEMA.md
#
# Including this file (as the test does) defines `schema_markdown()` + `SCHEMA_PATH` WITHOUT writing.

using PTCryspMC

const SCHEMA_PATH = normpath(joinpath(@__DIR__, "..", "docs", "SCHEMA.md"))

# One markdown table from a `(name, vector)` column list + its (unit, description) doc map. Type is
# the stored dtype (eltype of the column vector); a missing doc entry is a hard error.
function _column_table(io, columns, doc)
    println(io, "| Column | Type | Unit | Description |")
    println(io, "|--------|------|------|-------------|")
    for (name, v) in columns
        haskey(doc, name) || error("schema doc missing column '$name' — add it to the doc map")
        unit, desc = doc[name]
        println(io, "| `", name, "` | ", eltype(v), " | ", unit, " | ", desc, " |")
    end
end

"The full docs/SCHEMA.md as a string — deterministic, the single rendering used by both the writer below and the drift test."
function schema_markdown()::String
    io = IOBuffer()
    println(io, "# Output schema — PTCryspMC.jl")
    println(io)
    println(io, "**Generated from the code by `scripts/gen_schema.jl` — do not edit by hand.** ",
                "Regenerate with `julia --project=. scripts/gen_schema.jl`. Columns and types are ",
                "introspected from `singles_columns` / `coinc_columns` and the struct field types; ",
                "units and meaning come from the `singles_doc` / `coinc_doc` maps beside them. ",
                "`test/runtests.jl` fails if this file drifts from the code.")
    println(io)
    println(io, "Positions and energies are stored as quantized `Int16`: position = ",
                "`round(mm / $(XYZ_SCALE_MM))`, energy = `round(keV / $(E_SCALE_KEV))` — lossless at ",
                "detector resolution; decode with `decode_xyz` / `decode_e`. Times are `Float32` ",
                "nanoseconds **relative to the decay** (absolute = `event_time(activity, event)·1e9 + t`).")
    println(io)

    println(io, "## Singles — `prod/<tag>/singles.{h5,csv}`")
    println(io)
    println(io, "One row per detected photon (deposited in the ring; a miss writes nothing).")
    println(io)
    _column_table(io, singles_columns(SinglesBuffer()), PTCryspMC.singles_doc)
    println(io)

    println(io, "## LORs — `prod/<tag>/{lors_truth,lors_det,randoms}.h5`")
    println(io)
    println(io, "One row per coincidence (line of response). All three files share this schema: ",
                "`lors_truth` = same-annihilation pairs (truth, unsmeared); `lors_det` = the smeared, ",
                "energy/DT-selected detector list (truth ∪ randoms); `randoms` = cross-decay accidentals ",
                "(`truth = 2`, `dt_ns = NaN`).")
    println(io)
    _column_table(io, coinc_columns(CoincidenceBuffer()), PTCryspMC.coinc_doc)
    println(io)

    println(io, "### Truth code")
    println(io)
    println(io, "| Value | Meaning |")
    println(io, "|-------|---------|")
    println(io, "| ", Int(TRUTH_TRUE),    " | true |")
    println(io, "| ", Int(TRUTH_SCATTER), " | scatter |")
    println(io, "| ", Int(TRUTH_RANDOM),  " | random |")
    println(io)
    println(io, "A `truth = 1` (scatter) LOR is **single** scatter when `nscat1 + nscat2 == 1`, ",
                "**multiple** when `≥ 2`.")
    println(io)

    println(io, "## HDF5 root attributes (provenance)")
    println(io)
    println(io, "Set by the writers (not part of the column lists), carried for provenance and exact ",
                "regeneration of pruned singles. Representative keys:")
    println(io)
    println(io, "| Attribute | Files | Meaning |")
    println(io, "|-----------|-------|---------|")
    println(io, "| `scenario_tag` | all | run tag (config filename) |")
    println(io, "| `nrows` | all | row count |")
    println(io, "| `xyz_scale_mm`, `e_scale_keV` | all | quantization scales ($(XYZ_SCALE_MM) mm, $(E_SCALE_KEV) keV) |")
    println(io, "| `seed` | singles | transport seed; on LORs, the builder's smear/time seed |")
    println(io, "| `transport_seed`, `nchunks`, `nevents` | singles, lors_truth, lors_det | the recipe to regenerate the singles exactly |")
    println(io, "| `crystal`, `phantom_material` | all | materials |")
    println(io, "| `mode`, `has_randoms` | LORs | `truth`/`det`/`randoms`; whether randoms are merged |")
    println(io, "| `tau_ns` | lors_det, randoms | coincidence window |")
    println(io, "| `sigma_xyz_mm`, `eres`, `emin_keV`, `window_fwhm` | lors_det | detector response + energy cut |")
    println(io, "| `t_relative_to_decay`, `t0_s`, `t1_s`, `half_life_s`, `time_seed` | LORs | the activity clock for absolute time |")
    String(take!(io))
end

if abspath(PROGRAM_FILE) == @__FILE__
    write(SCHEMA_PATH, schema_markdown())
    println("wrote ", SCHEMA_PATH, " (", length(schema_markdown()), " bytes)")
end
