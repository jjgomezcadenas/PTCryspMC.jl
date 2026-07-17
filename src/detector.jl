# Detector response: how a truth hit is blurred into a measurement. The energy
# resolution is FWHM(E) = a·√(511 keV / E), with `a` the fractional FWHM at 511 keV
# (5% pure CsI, 7% CsI(Tl), 10% BGO); position is a Gaussian σ_xyz that includes the DOI.
# (Times come later, with the randoms pass.) See docs/PTCryspMC_phys.tex.

const FWHM_TO_SIGMA = 2.3548200450309493        # 2√(2 ln 2)

"""
    energy_fwhm(E, a) -> Float64

Absolute energy-resolution FWHM [keV] at deposited energy `E` [keV] for a detector with
fractional FWHM `a` at 511 keV: FWHM(E) = a·√(511·E) (so the *fractional* width is
a·√(511/E), and equals `a` at E = 511 keV).
"""
energy_fwhm(E::Real, a::Real)::Float64 = a * sqrt(511.0 * Float64(E))

"Gaussian σ [keV] of the energy measurement at energy `E`, resolution `a` (= FWHM/2.355)."
energy_sigma(E::Real, a::Real)::Float64 = energy_fwhm(E, a) / FWHM_TO_SIGMA

"""
    smear_energy(E, a, rng) -> Float64

Blur a deposited energy `E` [keV] by the detector energy resolution `a` (fractional FWHM
at 511 keV). `a ≤ 0` returns `E` unchanged (perfect resolution / truth).
"""
@inline smear_energy(E::Float64, a::Float64, rng::AbstractRNG)::Float64 =
    a > 0.0 ? E + energy_sigma(E, a) * randn(rng) : E

"""
    smear_position(p, σ, rng) -> NTuple{3,Float64}

Blur a 3-D position by an isotropic Gaussian of width `σ` (each coordinate independently).
`σ ≤ 0` returns `p` unchanged.
"""
@inline function smear_position(p, σ::Float64, rng::AbstractRNG)::NTuple{3,Float64}
    σ > 0.0 || return (Float64(p[1]), Float64(p[2]), Float64(p[3]))
    (p[1] + σ * randn(rng), p[2] + σ * randn(rng), p[3] + σ * randn(rng))
end

"""
    smear_position_mix(p, σ1, σ2, f, rng) -> NTuple{3,Float64}

Position model 2: blur each coordinate independently by a two-Gaussian core/tail mixture —
with probability `f` (the core fraction) a Gaussian of width `σ1`, else the tail `σ2`. The
core/tail draw is per coordinate (the model is defined by the per-coordinate marginal fits;
`data/pet_resolution_recipe.json`, averaged over x/y/z). `σ1 ≤ 0` returns `p` unchanged.
"""
@inline function smear_position_mix(p, σ1::Float64, σ2::Float64, f::Float64,
                                    rng::AbstractRNG)::NTuple{3,Float64}
    σ1 > 0.0 || return (Float64(p[1]), Float64(p[2]), Float64(p[3]))
    s(x) = x + (rand(rng) < f ? σ1 : σ2) * randn(rng)
    (s(p[1]), s(p[2]), s(p[3]))
end
