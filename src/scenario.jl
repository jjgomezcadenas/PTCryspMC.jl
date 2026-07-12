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
    _activity_distal_edge_mm(pools, weights; binw=2.0) -> Float64

The activity distal edge [mm, z]: the most-downstream (+z) point where the weighted activity
depth profile falls through half its peak. `weights[j]` scales isotope j's pooled points to the
true abundance (as in `write_activity_profile`), so the edge is the physical R50 of the total
activity. Used to centre the source on the range endpoint (`center_on="distal_edge"`).
"""
function _activity_distal_edge_mm(pools::Vector{Vector{NTuple{3,Float64}}},
                                  weights::Vector{Float64}; binw::Float64=2.0)::Float64
    zmin = Inf; zmax = -Inf
    for pool in pools, p in pool
        z = p[3] * 10
        z < zmin && (zmin = z); z > zmax && (zmax = z)
    end
    isfinite(zmin) || error("distal-edge centring: no emitter points")
    nb = max(2, ceil(Int, (zmax - zmin) / binw) + 1)
    h  = zeros(Float64, nb)
    for (j, pool) in enumerate(pools)
        w = weights[j]; w == 0.0 && continue
        for p in pool
            b = clamp(floor(Int, (p[3] * 10 - zmin) / binw) + 1, 1, nb)
            h[b] += w
        end
    end
    pk = argmax(h); half = h[pk] / 2
    edge = zmin + (pk - 0.5) * binw
    for b in pk:(nb - 1)
        if h[b] >= half && h[b + 1] < half
            z1   = zmin + (b - 0.5) * binw
            frac = (h[b] - half) / (h[b] - h[b + 1])
            edge = z1 + frac * binw
            break
        end
    end
    edge
end

"""
    load_scenario(dir, materials; budget="fast", dose_Gy=1.0, keep_escaped=false, center_on="") -> Scenario

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
                       keep_escaped::Bool=false, center_on::AbstractString="")::Scenario
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

    # Optional patient positioning: rigidly shift the source + phantom in z so a chosen reference
    # point sits at the ring centre (z=0). `center_on="distal_edge"` centres the activity distal
    # edge — the range endpoint — which maximises the acceptance of the LORs from that critical
    # region. The scanner is untouched; only the relative source/phantom placement moves, exactly
    # as a patient is positioned. Emitters and phantom shift together, so the inside/escaped
    # classification (done above, in the original frame) is preserved.
    z_offset_mm = 0.0
    if center_on == "distal_edge"
        weights = [length(pools[j]) > 0 ? n_exp[j] * f_inside[j] / length(pools[j]) : 0.0 for j in 1:n_iso]
        z_offset_mm = -_activity_distal_edge_mm(pools, weights)
    elseif !isempty(center_on)
        error("[source].center_on '$center_on' unknown (supported: \"distal_edge\", or omit)")
    end
    if z_offset_mm != 0.0
        off = z_offset_mm / 10                      # mm → cm
        for pool in pools, i in eachindex(pool)
            @inbounds pool[i] = (pool[i][1], pool[i][2], pool[i][3] + off)
        end
        p0 = phantom.position
        phantom = PhysicalVolume(phantom.logical, (p0[1], p0[2], p0[3] + off))
    end

    rmeta = _read_meta_row(joinpath(dir, "run_meta.csv"))
    sname = basename(rstrip(dir, '/'))
    prov = Dict{String,Any}(
        "scenario" => sname, "budget" => budget, "dose_Gy" => Float64(dose_Gy),
        "center_on" => center_on, "z_offset_mm" => z_offset_mm,
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
    write_activity_profile(path, scn, depth_dose_path)

Derive the binned true β⁺ activity(z) — the clean, detector-independent source curve the analysis
scores its reconstruction against — and write it to `path`. Bins are the exact z-frame of
`depth_dose.csv` (`z_mm` column), so activity(z) and dose(z) overlay directly; the outermost bin
edges run to ±∞ so every inside-phantom emitter is counted and each isotope's column integrates to
its expected decay count at the scenario's `dose_Gy`. Per isotope j the pool holds `N_kept_j` sampled
annihilation points (escaped positrons already dropped), representing the inside distribution from
which the shard materializes `M_j ~ Poisson(N_expected_j·f_inside_j)`; so each pooled point carries
weight `N_expected_j·f_inside_j / N_kept_j` — the same source scaling the LOR shard uses, hence the
profile composes with the LORs. Columns: `z_mm`, one per isotope (its `name`, in isotope-id order),
`total`. Writes a companion `<path stem>_meta.csv` (scenario/budget/dose/units).
"""
function write_activity_profile(path::AbstractString, scn::Scenario, depth_dose_path::AbstractString)
    dh, dr = _read_csv(depth_dose_path)
    zc = _col(dh, "z_mm")
    centers = [parse(Float64, r[zc]) for r in dr]     # mm, the dose z-frame
    nb = length(centers)
    nb >= 2 || error("depth_dose.csv needs ≥2 z bins (got $nb)")
    # interior edges (midpoints); the two outer bins run to ±∞ so nothing inside-phantom is lost.
    edges = [(centers[i] + centers[i + 1]) / 2 for i in 1:(nb - 1)]

    niso = length(scn.pools)
    acc  = [zeros(Float64, nb) for _ in 1:niso]        # acc[j][bin] = expected decays
    for j in 1:niso
        nkept = length(scn.pools[j])
        nkept == 0 && continue
        w = scn.n_expected[j] * scn.f_inside[j] / nkept
        col = acc[j]
        @inbounds for p in scn.pools[j]
            b = searchsortedlast(edges, p[3] * 10) + 1  # cm → mm; 1..nb
            col[b] += w
        end
    end

    names = [scn.isotopes[j].name for j in 1:niso]
    open(path, "w") do io
        println(io, "z_mm,", join(names, ","), ",total")
        for k in 1:nb
            tot = 0.0
            print(io, centers[k])
            for j in 1:niso
                v = acc[j][k]; tot += v
                print(io, ",", v)
            end
            println(io, ",", tot)
        end
    end
    # companion meta (self-describing): the scaling this profile carries.
    meta = replace(path, r"\.csv$" => "_meta.csv")
    open(meta, "w") do io
        println(io, "scenario,budget,dose_Gy,target_dose_Gy,t_meas_s,n_isotopes,z_frame,units,note")
        println(io, join((scn.name, scn.budget, scn.dose_Gy, scn.target_dose_Gy, scn.t_meas_s,
                          niso, "depth_dose.csv", "decays",
                          "expected(mean) counts at dose_Gy; dose-linear; escaped positrons excluded"), ","))
    end
    path
end

"""
    write_truth_bundle(scndir, truthdir, materials; budget="fast", dose_Gy=1.0, force=false)

Export the detector-independent scenario truth into `truthdir` (the `<scenario>/truth/` level of the
products tree, shared across scanners/crystals like `phantom/`). Five verbatim copies —
`depth_dose.csv` (→ dose-R80), `sobp_layers.csv` (+ `_meta`), `run_meta.csv`, and the budget's
`sampling_budget_<budget>.csv` (+ `_meta`) — plus the one derived `activity_profile_<budget>.csv`
(→ activity-R50, via `write_activity_profile`, which streams `emitters.csv`). Idempotent: existing
files are kept unless `force`. Called by `publish_prod` (once per scenario) and runnable stand-alone
to backfill an already-published master. Returns `truthdir`.
"""
function write_truth_bundle(scndir::AbstractString, truthdir::AbstractString,
                            materials::Dict{String,Material};
                            budget::AbstractString="fast", dose_Gy::Real=1.0, force::Bool=false)
    isdir(scndir) || error("scenario_dir not found: $scndir")
    mkpath(truthdir)
    copies = ["depth_dose.csv", "sobp_layers.csv", "sobp_layers_meta.csv", "run_meta.csv",
              "sampling_budget_$(budget).csv", "sampling_budget_$(budget)_meta.csv"]
    for f in copies
        src = joinpath(scndir, f)
        isfile(src) || error("truth source missing: $src")
        dst = joinpath(truthdir, f)
        (force || !isfile(dst)) && cp(src, dst; force=true)
    end
    prof = joinpath(truthdir, "activity_profile_$(budget).csv")
    if force || !isfile(prof)
        scn = load_scenario(scndir, materials; budget=budget, dose_Gy=dose_Gy)
        write_activity_profile(prof, scn, joinpath(scndir, "depth_dose.csv"))
    end
    truthdir
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
