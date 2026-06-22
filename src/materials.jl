# Materials: density + XCOM photon cross sections, and the macroscopic cross
# sections used in transport (photon-only). Because the XCOM table is per-compound,
# the macroscopic cross section is just Sigma = (mu/rho) * density [cm^-1] — no
# per-atom / element-mixing bookkeeping is needed here.

"""
    Material

A material's density and its XCOM photon cross sections (incoherent, photoelectric,
pair = nuclear + electron), with pre-logged grids for fast interpolation. Vacuum
(no `xcom`) has empty grids and zero cross sections.
"""
struct Material
    name::String
    density::Float64                  # g/cm^3
    # Scintillation / readout properties, keyed to the crystal (0 / empty for non-scintillators
    # like Water/Air/Vacuum). The crystal name pulls these — no extra config cards.
    light_yield::Float64              # scintillation photons per MeV
    scint_decay_ns::Vector{Float64}   # decay-component lifetimes [ns]
    scint_decay_w::Vector{Float64}    # component weights (sum to 1; one entry = single exponential)
    eres_a::Float64                   # energy resolution: fractional FWHM at 511 keV (0 = unset)
    pde::Float64                      # effective photon-detection efficiency for THIS crystal's
                                      # emission spectrum (per-crystal: depends on emission colour)
    E::Vector{Float64}                # nudged energy grid [MeV]
    log_E::Vector{Float64}
    incoherent::Vector{Float64}       # mu/rho [cm^2/g]
    photoelectric::Vector{Float64}
    pair::Vector{Float64}             # nuclear + electron
    log_incoherent::Vector{Float64}
    log_photoelectric::Vector{Float64}
    log_pair::Vector{Float64}
end

"Read and parse `materials.json` from `data_dir` into its raw Dict."
_read_materials_json(data_dir::AbstractString) =
    open(io -> JSON.parse(io), joinpath(data_dir, "materials.json"), "r")

"""
    _build_material(data_dir, name, d) -> Material

Build one `Material` from its `materials.json` entry `d`. Materials without an
`xcom` entry (e.g. Vacuum) get empty grids and zero cross sections.
"""
function _build_material(data_dir::AbstractString, name::AbstractString, d)::Material
    density = Float64(d["density"])
    yield   = Float64(get(d, "light_yield_per_MeV", 0.0))
    decay   = Float64.(get(d, "scint_decay_ns", Float64[]))
    w       = Float64.(get(d, "scint_decay_w", Float64[]))
    eres_a  = Float64(get(d, "eres_a", 0.0))
    pde     = Float64(get(d, "pde", 0.0))
    length(decay) == length(w) ||
        error("material '$name': scint_decay_ns and scint_decay_w must have equal length")
    (isempty(w) || isapprox(sum(w), 1.0; atol=1e-6)) ||
        error("material '$name': scint_decay_w must sum to 1 (got $(sum(w)))")

    xf = get(d, "xcom", nothing)
    if xf === nothing
        return Material(name, density, yield, decay, w, eres_a, pde, Float64[], Float64[],
                        Float64[], Float64[], Float64[], Float64[], Float64[], Float64[])
    end
    xc = load_xcom(joinpath(data_dir, xf))
    E = _prepare_xcom_energy(xc)
    pair = xc.pair_nuclear .+ xc.pair_electron
    Material(name, density, yield, decay, w, eres_a, pde,
             E, log.(E), xc.incoherent, xc.photoelectric, pair,
             prelog_data(xc.incoherent), prelog_data(xc.photoelectric), prelog_data(pair))
end

"""
    load_material(data_dir, name) -> Material

Load a single material by name from `materials.json` in `data_dir`, reading only
its XCOM table. Use this when a run needs one (or a few) named materials, so an
unrelated or unfinished entry in `materials.json` cannot break the run.
"""
function load_material(data_dir::AbstractString, name::AbstractString)::Material
    raw = _read_materials_json(data_dir)
    haskey(raw, name) || error("material '$name' not found in materials.json")
    _build_material(data_dir, name, raw[name])
end

"""
    load_materials(data_dir) -> Dict{String, Material}

Read `materials.json` from `data_dir` and build every material with an `xcom`
entry. Keys starting with `_` (documentation) are skipped.
"""
function load_materials(data_dir::AbstractString)::Dict{String,Material}
    raw = _read_materials_json(data_dir)
    mats = Dict{String,Material}()
    for (name, d) in raw
        startswith(name, "_") && continue
        mats[name] = _build_material(data_dir, name, d)
    end
    mats
end

"""
    sigma_macro(mat, E_MeV) -> (Sigma_C, Sigma_Ph, Sigma_P)

Macroscopic cross sections [cm^-1] for Compton (incoherent), photoelectric, and
pair: Sigma = (mu/rho)(E) * density. Pair is last as it is zero below threshold
(<= 511 keV). Vacuum returns zeros.
"""
function sigma_macro(mat::Material, E_MeV::Float64)::Tuple{Float64,Float64,Float64}
    isempty(mat.E) && return (0.0, 0.0, 0.0)
    lx = log(E_MeV)
    n = length(mat.E)
    lo = clamp(searchsortedlast(mat.E, E_MeV), 1, n - 1)
    ρ = mat.density
    ΣC  = interp_loglog_prelogged(lx, mat.log_E, mat.log_incoherent,    mat.incoherent,    lo) * ρ
    ΣPh = interp_loglog_prelogged(lx, mat.log_E, mat.log_photoelectric, mat.photoelectric, lo) * ρ
    ΣP  = interp_loglog_prelogged(lx, mat.log_E, mat.log_pair,          mat.pair,          lo) * ρ
    (ΣC, ΣPh, ΣP)
end

"Photon mean free path [cm] in `mat` at energy `E_MeV`."
function mfp(mat::Material, E_MeV::Float64)::Float64
    Σ = sum(sigma_macro(mat, E_MeV))
    Σ > 0.0 ? 1.0 / Σ : Inf
end
