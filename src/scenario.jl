# Reading a frozen ptcryspg4 scenario — the API (Proton Activity) source mode. Parses the
# handoff CSVs (phantom medium here; emitters + budget + isotopes with the scenario reader)
# into the engine's types. The scenario is the single source of truth: the phantom is built
# from its phantom_regions.csv, not a hand-written JSON. See dev/api_plan.md and
# docs/PTCryspMC_app.tex §3. All positions in the scenario are mm in the origin-centred,
# beam-+z world frame shared with emitters.csv; the engine works in cm, so lengths ÷ 10.

"""
    _read_csv(path) -> (header::Vector{String}, rows::Vector{Vector{String}})

Minimal reader for the small comma-separated scenario tables (unquoted). Returns the
header names and the data rows as trimmed strings; blank lines are skipped.
"""
function _read_csv(path::AbstractString)::Tuple{Vector{String},Vector{Vector{String}}}
    lines = filter(!isempty, strip.(readlines(path)))
    isempty(lines) && error("empty CSV: $path")
    header = String.(strip.(split(lines[1], ',')))
    rows   = [String.(strip.(split(l, ','))) for l in lines[2:end]]
    (header, rows)
end

"Index of column `name` in a `_read_csv` header, or a clear error listing what is present."
function _col(header::Vector{String}, name::AbstractString)::Int
    i = findfirst(==(name), header)
    i === nothing && error("column '$name' not found (have: $(join(header, ", ")))")
    i
end

"""
    load_phantom_regions(path, materials) -> PhysicalVolume

Build the phantom volume from a scenario's `phantom_regions.csv` (mm → cm). One row per
region: `solid` (`ellipsoid` | `cylinder`), the semi-axes/radius `a_mm,b_mm,c_mm` (cylinder:
`radius, radius, half-length`), the centre `cx_mm,cy_mm,cz_mm`, the `material` (keyed directly
into `materials`), and Euler angles. The region frame is the scenario world frame (origin-
centred, beam +z), co-registered with `emitters.csv`, so the centre becomes the placement
directly and the cylinder axis is z (the beam).

**Single-region only** for now: errors on >1 region (multi-region phantom deferred — the
navigator has one phantom leaf) and on any nonzero Euler angle (rotation deferred; the transform
belongs in `PhysicalVolume._to_local` when added). Errors clearly on an unknown material or an
unsupported solid rather than mis-loading.
"""
function load_phantom_regions(path::AbstractString, materials::Dict{String,Material})::PhysicalVolume
    header, rows = _read_csv(path)
    isempty(rows) && error("no regions in $path")
    length(rows) == 1 ||
        error("multi-region phantom deferred: $(length(rows)) regions in $path " *
              "(only single-region supported — the navigator has one phantom leaf)")
    r = rows[1]
    getf(name) = r[_col(header, name)]

    for e in ("euler_x_deg", "euler_y_deg", "euler_z_deg")
        parse(Float64, getf(e)) == 0.0 ||
            error("region rotation deferred: $e = $(getf(e)) (only axis-aligned regions supported)")
    end

    a = parse(Float64, getf("a_mm")) / 10   # mm → cm
    b = parse(Float64, getf("b_mm")) / 10
    c = parse(Float64, getf("c_mm")) / 10
    shape = getf("solid")
    sol = if shape == "ellipsoid"
        Ellipsoid(a, b, c)
    elseif shape == "cylinder"
        isapprox(a, b) ||
            error("cylinder region needs a_mm == b_mm (radius, radius); got $a, $b cm")
        Cylinder(a, c)                      # (radius, half-length); axis along z = beam
    else
        error("unsupported region solid '$shape' (ellipsoid, cylinder)")
    end

    matname = getf("material")
    haskey(materials, matname) ||
        error("region material '$matname' not in materials (add it to data/materials.json)")

    centre = (parse(Float64, getf("cx_mm")) / 10,
              parse(Float64, getf("cy_mm")) / 10,
              parse(Float64, getf("cz_mm")) / 10)
    PhysicalVolume(LogicalVolume(getf("region"), sol, materials[matname]), centre)
end
