# Emission source: where annihilations happen and how the back-to-back pair leaves.
#
# A Source draws an annihilation point and emits two ~511 keV photons. The point is
# drawn uniformly inside a phantom volume (source = attenuator: the same solid both
# emits and attenuates — see docs/PTCryspMC_app.tex). The two photons are back to
# back up to a small acollinearity (~0.5 deg FWHM): they are not exactly 180 deg apart.

"An isotropic unit direction (uniform on the sphere)."
@inline function rand_direction(rng::AbstractRNG)
    c = 2.0 * rand(rng) - 1.0
    ϕ = 2π * rand(rng)
    s = sqrt(max(0.0, 1.0 - c^2))
    (s * cos(ϕ), s * sin(ϕ), c)
end

"""
    sample_point_in(solid, rng) -> NTuple{3,Float64}

A point drawn uniformly by volume inside `solid`, in the solid's local frame.
"""
function sample_point_in(c::Cylinder, rng::AbstractRNG)
    r = c.radius_cm * sqrt(rand(rng))                 # uniform in the disc
    ϕ = 2π * rand(rng)
    z = (2.0 * rand(rng) - 1.0) * c.half_length_cm
    (r * cos(ϕ), r * sin(ϕ), z)
end

function sample_point_in(s::Sphere, rng::AbstractRNG)
    r = s.radius_cm * cbrt(rand(rng))                 # uniform in the ball
    u = rand_direction(rng)
    (r * u[1], r * u[2], r * u[3])
end

function sample_point_in(b::Box, rng::AbstractRNG)
    ((2.0 * rand(rng) - 1.0) * b.half_x_cm,
     (2.0 * rand(rng) - 1.0) * b.half_y_cm,
     (2.0 * rand(rng) - 1.0) * b.half_z_cm)
end

abstract type Source end

"""
    UniformVolumeSource(volume)

A source filling a placed phantom `volume` (a `PhysicalVolume`) with uniform activity:
annihilation points are drawn uniformly inside its solid and placed in the world frame.
The same volume is also the attenuator (its material), so a Vacuum/Air phantom emits
without attenuating (the non-attenuated reference) and a Water phantom attenuates.
"""
struct UniformVolumeSource{S<:Solid} <: Source
    volume::PhysicalVolume{S}
end

"A world-frame annihilation point drawn from the source."
function sample_position(src::UniformVolumeSource, rng::AbstractRNG)
    p = sample_point_in(solid(src.volume), rng)
    o = src.volume.position
    (p[1] + o[1], p[2] + o[2], p[3] + o[3])
end

"""
    PointSource(position)

A point source at a fixed world-frame `position` (e.g. the phantom centre) — the
degenerate phantom, useful for quick checks. Goes through the same `emit_pair`
acollinearity machinery as a volume source.
"""
struct PointSource <: Source
    position::NTuple{3,Float64}
end

sample_position(src::PointSource, ::AbstractRNG) = src.position

# Tilt the unit vector `axis` by a small 2-D Gaussian deviation in its transverse plane:
# each transverse component ~ N(0, σ), so each projected angle has the given FWHM. As
# σ -> 0 the result is `axis` exactly. Used for the back-to-back acollinearity.
@inline function _acollinear(axis, σ::Float64, rng::AbstractRNG)
    σ <= 0.0 && return axis
    ax = σ * randn(rng); ay = σ * randn(rng)
    s2 = ax^2 + ay^2
    # Unit local direction with transverse components (ax, ay). For the normal small-angle
    # case use sqrt(1-s2) as the z-cosine; fall back to a renormalised tilt if σ is so
    # large that s2 ≥ 1 (never happens for sub-degree acollinearity, but keep it safe).
    local_dir = s2 < 1.0 ? (ax, ay, sqrt(1.0 - s2)) : (ax, ay, 1.0) ./ sqrt(s2 + 1.0)
    # The tuple twin returns an NTuple directly — no per-annihilation heap Vector (this runs once
    # per decay on the production MT loop, so the allocation would otherwise defeat the alloc-free
    # transport's GC-quiet goal).
    rotate_to_global_t(local_dir[1], local_dir[2], local_dir[3], axis)
end

"""
    emit_pair(src, rng; acol_fwhm_deg=0.5) -> (pos, dir1, dir2)

One annihilation from `src`: a world-frame point and two ~back-to-back unit directions.
`dir1` is isotropic; `dir2` is the partner of `-dir1` tilted by a small Gaussian
acollinearity of `acol_fwhm_deg` FWHM per transverse projection (0 = exactly 180 deg).
"""
function emit_pair(src::Source, rng::AbstractRNG; acol_fwhm_deg::Float64=0.5)
    pos = sample_position(src, rng)
    d1  = rand_direction(rng)
    σ   = deg2rad(acol_fwhm_deg) / 2.3548200450309493     # FWHM -> sigma (2*sqrt(2 ln 2))
    d2  = _acollinear((-d1[1], -d1[2], -d1[3]), σ, rng)
    (pos, d1, d2)
end
