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

"Read a one-row CSV (run_meta, a budget meta) into a `String => String` field map."
function _read_meta_row(path::AbstractString)::Dict{String,String}
    header, rows = _read_csv(path)
    isempty(rows) && error("no data row in $path")
    Dict(header[i] => rows[1][i] for i in eachindex(header))
end

# Stream emitters.csv (the big file — millions of rows) into per-isotope pools of annihilation
# points (`anh_*`, mm → cm), one push per row so the whole file is never held as strings. Points
# outside the phantom are the escaped/surface-pinned positrons (they annihilate far away in air,
# not at the pinned point) — DROPPED unless `keep_escaped`. Returns the pools and the per-isotope
# total (kept + dropped), so the caller can form the inside fraction.
function _read_emitter_pools(path::AbstractString, n_iso::Int, phantom, keep_escaped::Bool)
    pools   = [NTuple{3,Float64}[] for _ in 1:n_iso]
    n_total = zeros(Int, n_iso)
    open(path) do io
        header = String.(strip.(split(strip(readline(io)), ',')))
        iid_c = _col(header, "isotope_id")
        ax_c  = _col(header, "anh_x_mm"); ay_c = _col(header, "anh_y_mm"); az_c = _col(header, "anh_z_mm")
        for line in eachline(io)
            isempty(line) && continue
            f = split(line, ',')
            iid = parse(Int, f[iid_c])
            (0 <= iid < n_iso) || error("emitter isotope_id $iid out of range [0, $(n_iso-1)] in $path")
            p = (parse(Float64, f[ax_c]) / 10, parse(Float64, f[ay_c]) / 10, parse(Float64, f[az_c]) / 10)
            n_total[iid + 1] += 1
            (keep_escaped || is_inside(phantom, p)) && push!(pools[iid + 1], p)
        end
    end
    (pools, n_total)
end

"""
    Scenario

A parsed frozen `ptcryspg4` scenario (the API source). Holds the phantom (the single source of
truth, from `phantom_regions.csv`), the per-isotope annihilation-point pools (cm, escaped positrons
dropped), the isotopes, the per-isotope expected decay budget `N_expected` (rescaled to `dose_Gy`),
the measurement window `t_meas_s` (the truncated-exponential timing window), the per-isotope kept
fraction `f_inside` and dropped count (the escaped positrons), and the provenance to stamp into
outputs. Pools/isotopes/budget/f_inside are indexed by `isotope_id + 1`.
"""
struct Scenario
    name::String
    budget::String
    phantom::PhysicalVolume
    pools::Vector{Vector{NTuple{3,Float64}}}
    isotopes::Vector{Isotope}
    n_expected::Vector{Float64}
    t_meas_s::Float64
    dose_Gy::Float64
    target_dose_Gy::Float64
    f_inside::Vector{Float64}
    n_dropped::Vector{Int}
    provenance::Dict{String,Any}
end

"""
    load_scenario(dir, materials; budget="fast", dose_Gy=1.0, keep_escaped=false) -> Scenario

Read a frozen `ptcryspg4` scenario directory into a `Scenario`. Builds the phantom from the
scenario's own `phantom_regions.csv`; streams `emitters.csv` into per-isotope `anh` pools (mm → cm),
dropping escaped positrons (outside the phantom) unless `keep_escaped`; reads `isotopes.csv`,
`sampling_budget_<budget>.csv` (+ meta) and `run_meta.csv`. `N_expected` is rescaled linearly from
the budget's reference dose to `dose_Gy`. The escaped loss is *reported* here (`f_inside`,
`n_dropped`); the source (`APISource`, step 5) applies it as `M_j ~ Poisson(N_expected_j · f_inside_j)`.
`isotopes.csv` carries no β⁺ (the budget already counts annihilations), so `Isotope.beta_plus = 1`
here and is unused in API mode.
"""
function load_scenario(dir::AbstractString, materials::Dict{String,Material};
                       budget::AbstractString="fast", dose_Gy::Real=1.0,
                       keep_escaped::Bool=false)::Scenario
    phantom = load_phantom_regions(joinpath(dir, "phantom_regions.csv"), materials)

    ih, ir = _read_csv(joinpath(dir, "isotopes.csv"))
    idc = _col(ih, "isotope_id"); nc = _col(ih, "name"); hc = _col(ih, "half_life_s")
    n_iso = maximum(parse(Int, r[idc]) for r in ir) + 1
    isos  = Vector{Isotope}(undef, n_iso)
    for r in ir
        id = parse(Int, r[idc])
        isos[id + 1] = Isotope(r[nc], parse(Float64, r[hc]), 1.0)   # β⁺ unused in API
    end

    pools, n_total = _read_emitter_pools(joinpath(dir, "emitters.csv"), n_iso, phantom, keep_escaped)
    f_inside  = [n_total[j] > 0 ? length(pools[j]) / n_total[j] : 1.0 for j in 1:n_iso]
    n_dropped = [n_total[j] - length(pools[j]) for j in 1:n_iso]

    bh, br = _read_csv(joinpath(dir, "sampling_budget_$(budget).csv"))
    bidc = _col(bh, "isotope_id"); nec = _col(bh, "N_expected")
    bmeta   = _read_meta_row(joinpath(dir, "sampling_budget_$(budget)_meta.csv"))
    ref_dose = parse(Float64, bmeta["dose_Gy"])
    t_meas   = parse(Float64, bmeta["t_meas_s"])
    tgt_dose = parse(Float64, bmeta["target_dose_Gy"])
    scale = Float64(dose_Gy) / ref_dose
    n_exp = zeros(Float64, n_iso)
    for r in br
        n_exp[parse(Int, r[bidc]) + 1] = parse(Float64, r[nec]) * scale
    end

    rmeta = _read_meta_row(joinpath(dir, "run_meta.csv"))
    sname = basename(rstrip(dir, '/'))
    prov = Dict{String,Any}(
        "scenario" => sname, "budget" => budget, "dose_Gy" => Float64(dose_Gy),
        "geometry" => get(rmeta, "geometry", ""), "phantom_material" => get(rmeta, "phantom_material", ""),
        "geant4_version" => get(rmeta, "geant4_version", ""), "physics_list" => get(rmeta, "physics_list", ""),
        "upstream_seed" => get(rmeta, "random_seed", ""), "n_protons" => get(rmeta, "n_protons", ""),
        "keep_escaped" => keep_escaped, "prompt_gamma_modeled" => false,
        "n_escaped_dropped" => sum(n_dropped))

    Scenario(sname, budget, phantom, pools, isos, n_exp, t_meas,
             Float64(dose_Gy), tgt_dose, f_inside, n_dropped, prov)
end

# Poisson draw without a Distributions dependency: Knuth for a small mean (exact, O(mean), no
# exp(-mean) underflow below ~700), the normal approximation for a large mean (≥ 50 → the error is
# negligible and Knuth would be needlessly long). Only a handful of draws per realization (one per
# isotope), so correctness across the mean range — not speed — is the point.
function _rand_poisson(rng::AbstractRNG, mean::Float64)::Int
    mean <= 0.0 && return 0
    if mean < 50.0
        L = exp(-mean); k = 0; p = 1.0
        while true
            k += 1; p *= rand(rng)
            p <= L && return k - 1
        end
    else
        n = round(Int, mean + sqrt(mean) * randn(rng))
        return n < 0 ? 0 : n
    end
end

"""
    materialize_api_source(scn; master_seed=1, realization=0) -> APISource

Phase 1 of the API source: draw the whole annihilation-point array from a `Scenario`, seeded ONLY
by `(master_seed, realization)` — independent of the transport chunking, so every detector config
sees the identical source. For each isotope j draw `M_j ~ Poisson(N_expected_j · f_inside_j)` (the
escaped-positron loss folded in), then sample `M_j` points with replacement from that isotope's pool.
Events are numbered `1..N = 1..ΣM_j` in isotope-id order.
"""
function materialize_api_source(scn::Scenario; master_seed::Integer=1, realization::Integer=0)::APISource
    rng  = MersenneTwister(UInt64(master_seed) + UInt64(realization))
    niso = length(scn.pools)
    M = Vector{Int}(undef, niso)
    for j in 1:niso
        M[j] = _rand_poisson(rng, scn.n_expected[j] * scn.f_inside[j])
        (M[j] == 0 || !isempty(scn.pools[j])) ||
            error("isotope $(scn.isotopes[j].name): M_j=$(M[j]) but empty pool (cannot sample)")
    end
    N = sum(M)
    points  = Vector{NTuple{3,Float64}}(undef, N)
    isotope = Vector{Int8}(undef, N)
    e = 0
    for j in 1:niso
        pool = scn.pools[j]; npool = length(pool)
        @inbounds for _ in 1:M[j]
            e += 1
            points[e]  = pool[rand(rng, 1:npool)]
            isotope[e] = Int8(j - 1)
        end
    end
    lambdas = [log(2.0) / iso.half_life_s for iso in scn.isotopes]
    APISource(points, isotope, lambdas)
end

"""
    scenario_activity_models(scn; seed=1234) -> Vector{ActivityModel}

Per-isotope activity models for the API randoms timing: one `ActivityModel` per isotope, sharing
the window [0, t_meas] and `seed`, differing only in the decay constant λ (from each isotope's
half-life). Production ends before the measurement window opens, so the decay time within the
window is a pure per-isotope truncated exponential on [0, t_meas] — exactly this model. Index by
`isotope_id + 1`; the randoms pass picks the model by each single's isotope column.
"""
function scenario_activity_models(scn::Scenario; seed::Integer=1234)::Vector{ActivityModel}
    [ActivityModel(; t0=0.0, t1=scn.t_meas_s, half_life_s=iso.half_life_s, seed=seed)
     for iso in scn.isotopes]
end
