# Materials: density + XCOM photon cross sections, and the macroscopic cross
# sections used in transport. Adapted from LXeMC (src/materials.jl), trimmed to
# photon-only (no ESTAR / bremsstrahlung). Because the XCOM table is per-compound,
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
    E::Vector{Float64}                # nudged energy grid [MeV]
    log_E::Vector{Float64}
    incoherent::Vector{Float64}       # mu/rho [cm^2/g]
    photoelectric::Vector{Float64}
    pair::Vector{Float64}             # nuclear + electron
    log_incoherent::Vector{Float64}
    log_photoelectric::Vector{Float64}
    log_pair::Vector{Float64}
end

"""
    load_materials(data_dir) -> Dict{String, Material}

Read `materials.json` from `data_dir`; for each material with an `xcom` entry,
load its NIST XCOM table.
"""
function load_materials(data_dir::AbstractString)::Dict{String,Material}
    raw = open(joinpath(data_dir, "materials.json"), "r") do io
        JSON.parse(io)
    end
    mats = Dict{String,Material}()
    for (name, d) in raw
        startswith(name, "_") && continue
        density = Float64(d["density"])
        xf = get(d, "xcom", nothing)
        if xf === nothing
            mats[name] = Material(name, density, Float64[], Float64[],
                                  Float64[], Float64[], Float64[],
                                  Float64[], Float64[], Float64[])
            continue
        end
        xc = load_xcom(joinpath(data_dir, xf))
        E = _prepare_xcom_energy(xc)
        pair = xc.pair_nuclear .+ xc.pair_electron
        mats[name] = Material(name, density, E, log.(E),
                              xc.incoherent, xc.photoelectric, pair,
                              prelog_data(xc.incoherent),
                              prelog_data(xc.photoelectric),
                              prelog_data(pair))
    end
    mats
end

"""
    sigma_macro(mat, E_MeV) -> (Sigma_C, Sigma_P, Sigma_Ph)

Macroscopic cross sections [cm^-1] for Compton (incoherent), pair, and
photoelectric: Sigma = (mu/rho)(E) * density. Vacuum returns zeros.
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
    (ΣC, ΣP, ΣPh)
end

"Photon mean free path [cm] in `mat` at energy `E_MeV`."
function mfp(mat::Material, E_MeV::Float64)::Float64
    Σ = sum(sigma_macro(mat, E_MeV))
    Σ > 0.0 ? 1.0 / Σ : Inf
end
