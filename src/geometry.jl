# Geometry: a Geant4-style hierarchy with three kinds of object.
#
#   Solid          — pure shape, centred at the origin (local frame). No material,
#                    no placement. Implements the ray interface (is_inside,
#                    distance_to_entry, distance_to_exit, volume).
#   LogicalVolume  — a solid + the material it is made of. Reusable; no placement.
#   PhysicalVolume — a logical volume placed at a position (the world frame). It
#                    transforms a world ray into the solid's local frame once, then
#                    delegates to the solid.
#
# Geant4 semantics: material lives on the logical volume, placement on the physical
# volume, and the world->local transform is factored into the physical volume rather
# than copied into every solid.

# Surface tolerances. SURFACE_EPS [cm] (1 pm) rejects the start surface so a ray
# does not immediately re-hit the boundary it sits on; PARALLEL_EPS guards the
# divide-by-zero when a ray is parallel to a surface (axial vs. lateral).
const SURFACE_EPS  = 1e-10
const PARALLEL_EPS = 1e-20

# =====================================================================
# Solids (pure shape, local frame centred at the origin)
# =====================================================================

"""
    Solid

Abstract supertype for geometric shapes, defined in their own local frame
(centred at the origin). Every concrete `Solid` implements, in that frame:

  - `is_inside(s, p) -> Bool`
  - `distance_to_exit(pos, dir, s) -> Float64`
  - `distance_to_entry(pos, dir, s) -> Float64`
  - `volume(s) -> Float64`  (cm^3)
"""
abstract type Solid end

"A cylinder, axis along z, centred at the origin [cm]."
struct Cylinder <: Solid
    radius_cm::Float64
    half_length_cm::Float64
end

volume(c::Cylinder)::Float64 = π * c.radius_cm^2 * 2.0 * c.half_length_cm

function is_inside(c::Cylinder, p)::Bool
    (p[1]^2 + p[2]^2 <= c.radius_cm^2) && (abs(p[3]) <= c.half_length_cm)
end

"""
    _surface_crossings(pos, dir, c) -> (t_near, t_far)

Nearest and farthest distances [cm] at which the ray from `pos` along `dir`
crosses the cylinder's surface (lateral wall or either cap), in the local frame
and counting only crossings ahead (t > SURFACE_EPS). `t_near = Inf` /
`t_far = -Inf` when there is no such crossing. Boundary tests are inclusive
(`<=`) so a ray through the exact rim is claimed, never dropped. Allocation-free:
this is the geometry hot path called once per transport step.
"""
@inline function _surface_crossings(pos, dir, c::Cylinder)::Tuple{Float64,Float64}
    R = c.radius_cm; H = c.half_length_cm
    dx = pos[1]; dy = pos[2]; dz = pos[3]
    t_near = Inf; t_far = -Inf

    # Lateral wall: |(d + t·dir)_xy| = R, a quadratic in t.
    a = dir[1]^2 + dir[2]^2
    if a > PARALLEL_EPS
        b = 2.0 * (dx * dir[1] + dy * dir[2])
        cc = dx^2 + dy^2 - R^2
        disc = b^2 - 4.0 * a * cc
        if disc >= 0.0
            sq = sqrt(disc)
            for t in ((-b + sq) / (2a), (-b - sq) / (2a))
                if t > SURFACE_EPS && abs(dz + t * dir[3]) <= H
                    t_near = min(t_near, t); t_far = max(t_far, t)
                end
            end
        end
    end

    # End caps: z = ±H planes, accepted where the crossing falls within radius R.
    if abs(dir[3]) > PARALLEL_EPS
        for zf in (H, -H)
            t = (zf - dz) / dir[3]
            if t > SURFACE_EPS
                rx = dx + t * dir[1]; ry = dy + t * dir[2]
                (rx^2 + ry^2 <= R^2) && (t_near = min(t_near, t); t_far = max(t_far, t))
            end
        end
    end
    (t_near, t_far)
end

"""
    distance_to_exit(pos, dir, s) -> Float64

Distance [cm] from `pos` along `dir` to where the ray leaves the solid — the
farthest forward surface crossing. For an interior `pos` this is the single exit
point; `Inf` if the ray never meets the solid.
"""
function distance_to_exit(pos, dir, c::Cylinder)::Float64
    _, t_far = _surface_crossings(pos, dir, c)
    t_far > 0.0 ? t_far : Inf      # t_far is already pre-filtered to > SURFACE_EPS in the helper
end

"""
    distance_to_entry(pos, dir, s) -> Float64

Distance [cm] from an exterior `pos` along `dir` to where the ray first enters
the solid. `Inf` if the ray misses, or if `pos` is already inside (one forward
crossing only — that is an exit, not an entry).
"""
function distance_to_entry(pos, dir, c::Cylinder)::Float64
    t_near, t_far = _surface_crossings(pos, dir, c)
    t_near < t_far ? t_near : Inf
end

"An axis-aligned box, half-widths along x, y, z [cm], centred at the origin."
struct Box <: Solid
    half_x_cm::Float64
    half_y_cm::Float64
    half_z_cm::Float64
end

volume(b::Box)::Float64 = 8.0 * b.half_x_cm * b.half_y_cm * b.half_z_cm

function is_inside(b::Box, p)::Bool
    (abs(p[1]) <= b.half_x_cm) && (abs(p[2]) <= b.half_y_cm) && (abs(p[3]) <= b.half_z_cm)
end

"""
    _slab_crossings(pos, dir, b) -> (t_near, t_far)

Ray-box intersection by the slab method, local frame: the entry and exit distances
[cm] of the ray through the three axis-aligned slabs. `(Inf, -Inf)` when the ray
misses (or runs parallel to and outside a slab). Allocation-free hot path.
"""
@inline function _slab_crossings(pos, dir, b::Box)::Tuple{Float64,Float64}
    t_near = -Inf; t_far = Inf
    @inbounds for i in 1:3
        h = i == 1 ? b.half_x_cm : i == 2 ? b.half_y_cm : b.half_z_cm
        p = pos[i]; d = dir[i]
        if abs(d) > PARALLEL_EPS
            t1 = (-h - p) / d; t2 = (h - p) / d
            lo = min(t1, t2); hi = max(t1, t2)
            t_near = max(t_near, lo); t_far = min(t_far, hi)
        elseif p < -h || p > h
            return (Inf, -Inf)        # parallel to this slab and outside it: no hit
        end
    end
    t_far < t_near ? (Inf, -Inf) : (t_near, t_far)
end

function distance_to_exit(pos, dir, b::Box)::Float64
    _, t_far = _slab_crossings(pos, dir, b)
    t_far > SURFACE_EPS ? t_far : Inf
end

function distance_to_entry(pos, dir, b::Box)::Float64
    t_near, t_far = _slab_crossings(pos, dir, b)
    (t_near > SURFACE_EPS && t_near < t_far) ? t_near : Inf
end

"A sphere, centred at the origin [cm]."
struct Sphere <: Solid
    radius_cm::Float64
end

volume(s::Sphere)::Float64 = (4.0 / 3.0) * π * s.radius_cm^3

is_inside(s::Sphere, p)::Bool = p[1]^2 + p[2]^2 + p[3]^2 <= s.radius_cm^2

"""
    _sphere_crossings(pos, dir, s) -> (t_near, t_far)

Nearest and farthest distances [cm] at which the ray from `pos` along `dir` crosses
the sphere's surface (the two roots of ‖pos + t·dir‖ = R), in the local frame and
counting only crossings ahead (t > SURFACE_EPS). `(Inf, -Inf)` when there is none.
Allocation-free hot path.
"""
@inline function _sphere_crossings(pos, dir, s::Sphere)::Tuple{Float64,Float64}
    R = s.radius_cm
    a  = dir[1]^2 + dir[2]^2 + dir[3]^2
    t_near = Inf; t_far = -Inf
    a > PARALLEL_EPS || return (t_near, t_far)
    b  = 2.0 * (pos[1] * dir[1] + pos[2] * dir[2] + pos[3] * dir[3])
    cc = pos[1]^2 + pos[2]^2 + pos[3]^2 - R^2
    disc = b^2 - 4.0 * a * cc
    if disc >= 0.0
        sq = sqrt(disc)
        for t in ((-b + sq) / (2a), (-b - sq) / (2a))
            if t > SURFACE_EPS
                t_near = min(t_near, t); t_far = max(t_far, t)
            end
        end
    end
    (t_near, t_far)
end

function distance_to_exit(pos, dir, s::Sphere)::Float64
    _, t_far = _sphere_crossings(pos, dir, s)
    t_far > 0.0 ? t_far : Inf
end

function distance_to_entry(pos, dir, s::Sphere)::Float64
    t_near, t_far = _sphere_crossings(pos, dir, s)
    t_near < t_far ? t_near : Inf
end

"A cylindrical shell (hollow cylinder), axis along z, centred at the origin [cm]."
struct CylShell <: Solid
    r_inner_cm::Float64
    wall_thickness_cm::Float64
    half_length_cm::Float64
end

r_outer(cs::CylShell)::Float64 = cs.r_inner_cm + cs.wall_thickness_cm
volume(cs::CylShell)::Float64 =
    π * (r_outer(cs)^2 - cs.r_inner_cm^2) * 2.0 * cs.half_length_cm

function is_inside(cs::CylShell, p)::Bool
    r2 = p[1]^2 + p[2]^2
    (cs.r_inner_cm^2 <= r2 <= r_outer(cs)^2) && (abs(p[3]) <= cs.half_length_cm)
end

# Parameter interval [lo, hi] where the ray is inside the disc ρ² ≤ R² (the
# transverse quadratic a t² + b t + (pperp² − R²) = 0). Empty → (Inf, -Inf).
@inline function _disc_interval(a::Float64, b::Float64, pperp2::Float64, R2::Float64)
    if a > PARALLEL_EPS
        c = pperp2 - R2
        Δ = b * b - 4.0 * a * c
        Δ < 0.0 && return (Inf, -Inf)
        s = sqrt(Δ)
        return ((-b - s) / (2a), (-b + s) / (2a))
    end
    pperp2 <= R2 ? (-Inf, Inf) : (Inf, -Inf)   # parallel to axis: ρ constant
end

# Parameter interval [lo, hi] where the ray is inside the slab |z| ≤ H.
@inline function _z_interval(pz::Float64, dz::Float64, H::Float64)
    if abs(dz) > PARALLEL_EPS
        t1 = (-H - pz) / dz; t2 = (H - pz) / dz
        return (min(t1, t2), max(t1, t2))
    end
    abs(pz) <= H ? (-Inf, Inf) : (Inf, -Inf)   # parallel to caps
end

"""
    _shell_intervals(pos, dir, cs) -> (s1, e1, s2, e2)

The up-to-two parameter intervals [s1,e1], [s2,e2] where the ray lies in the shell
wall, S = (outer ∩ z-slab) ∖ bore. Absent pieces are returned as (Inf, -Inf). See
docs/navigation.tex. Allocation-free hot path.
"""
@inline function _shell_intervals(pos, dir, cs::CylShell)
    a = dir[1]^2 + dir[2]^2
    b = 2.0 * (pos[1] * dir[1] + pos[2] * dir[2])
    pperp2 = pos[1]^2 + pos[2]^2

    olo, ohi = _disc_interval(a, b, pperp2, r_outer(cs)^2)   # inside outer disc
    zlo, zhi = _z_interval(pos[3], dir[3], cs.half_length_cm)
    ta = max(olo, zlo); tb = min(ohi, zhi)                   # outer ∩ slab
    ta > tb && return (Inf, -Inf, Inf, -Inf)                 # never in the wall

    blo, bhi = _disc_interval(a, b, pperp2, cs.r_inner_cm^2) # inside the bore
    if blo > bhi || bhi <= ta || blo >= tb
        return (ta, tb, Inf, -Inf)                           # bore irrelevant: one piece
    end
    # bore (blo, bhi) carves [ta, tb] into the wall before it and after it
    s1, e1 = blo <= ta ? (Inf, -Inf) : (ta,  blo)
    s2, e2 = bhi >= tb ? (Inf, -Inf) : (bhi, tb)
    (s1, e1, s2, e2)
end

"""
    distance_to_exit(pos, dir, cs::CylShell) -> Float64

Distance to where the ray leaves the wall — the end of the shell interval the
photon is in (outer wall, inner wall into the bore, or a cap). `Inf` if `pos` is
not in the wall.
"""
function distance_to_exit(pos, dir, cs::CylShell)::Float64
    s1, e1, s2, e2 = _shell_intervals(pos, dir, cs)
    (s1 <= SURFACE_EPS && e1 > SURFACE_EPS) && return e1
    (s2 <= SURFACE_EPS && e2 > SURFACE_EPS) && return e2
    Inf
end

"""
    distance_to_entry(pos, dir, cs::CylShell) -> Float64

Distance to where the ray first enters the wall — the nearest shell-interval start
ahead. `Inf` if the ray misses the wall (or only grazes the bore).
"""
function distance_to_entry(pos, dir, cs::CylShell)::Float64
    s1, e1, s2, e2 = _shell_intervals(pos, dir, cs)
    t = Inf
    s1 > SURFACE_EPS && (t = min(t, s1))
    s2 > SURFACE_EPS && (t = min(t, s2))
    t
end

# =====================================================================
# Logical volumes (solid + material)
# =====================================================================

"A solid together with the material it is made of. No placement."
struct LogicalVolume{S<:Solid}
    name::String
    solid::S
    material::Material
end

name(lv::LogicalVolume)           = lv.name
solid(lv::LogicalVolume)          = lv.solid
material(lv::LogicalVolume)       = lv.material
volume(lv::LogicalVolume)::Float64 = volume(lv.solid)
"Mass [g] of the logical volume = density · volume."
mass(lv::LogicalVolume)::Float64 = lv.material.density * volume(lv.solid)

# =====================================================================
# Physical volumes (logical volume placed in the world)
# =====================================================================

"A logical volume placed at `position` in the world frame [cm]."
struct PhysicalVolume{S<:Solid}
    logical::LogicalVolume{S}
    position::NTuple{3,Float64}
end

solid(pv::PhysicalVolume)        = pv.logical.solid
material(pv::PhysicalVolume)     = pv.logical.material
name(pv::PhysicalVolume)         = pv.logical.name
volume(pv::PhysicalVolume)       = volume(pv.logical)
mass(pv::PhysicalVolume)         = mass(pv.logical)

# World -> local: subtract the placement, then delegate to the solid. (Rotation,
# when added, also rotates `dir` here — the single place it belongs.)
@inline _to_local(pv::PhysicalVolume, p) =
    (p[1] - pv.position[1], p[2] - pv.position[2], p[3] - pv.position[3])

is_inside(pv::PhysicalVolume, p)::Bool = is_inside(solid(pv), _to_local(pv, p))
distance_to_exit(pos, dir, pv::PhysicalVolume)::Float64 =
    distance_to_exit(_to_local(pv, pos), dir, solid(pv))
distance_to_entry(pos, dir, pv::PhysicalVolume)::Float64 =
    distance_to_entry(_to_local(pv, pos), dir, solid(pv))

# =====================================================================
# Detector scanner: a placed CylShell + its (φ, z) block / wheel grid
# =====================================================================

"""
    Scanner

The detector ring: a placed `CylShell` physical volume plus its readout
segmentation into `n_phi` azimuthal blocks × `n_z` axial wheels. A block is one
(φ, z) cell spanning the full wall depth (the radial depth is the DOI). See
docs/navigation.tex.
"""
struct Scanner
    volume::PhysicalVolume{CylShell}   # Deferred (dev/status.md #2): name collides with volume(); rename → pvol/shell
    n_phi::Int
    n_z::Int
end

"Number of readout blocks (= n_phi · n_z)."
nblocks(sc::Scanner)::Int = sc.n_phi * sc.n_z

"""
    block_index(sc, p) -> (iφ, iz)

The azimuthal block and axial wheel indices of a world-frame point `p` in the
scanner shell: φ = atan2(y,x) sectored into n_phi, z sectored into n_z. φ = 0 on
the +x axis, increasing counter-clockwise; the top edges (φ = 2π, z = +H) clamp
to the last cell.
"""
function block_index(sc::Scanner, p)::Tuple{Int,Int}
    lp = _to_local(sc.volume, p)
    H  = solid(sc.volume).half_length_cm
    φ  = mod(atan(lp[2], lp[1]), 2π)
    iφ = clamp(floor(Int, φ / (2π / sc.n_phi)), 0, sc.n_phi - 1)
    iz = clamp(floor(Int, (lp[3] + H) / (2H / sc.n_z)), 0, sc.n_z - 1)
    (iφ, iz)
end

"Linear block id of a world-frame point: iz·n_phi + iφ, in 0 … nblocks−1."
function block_id(sc::Scanner, p)::Int
    iφ, iz = block_index(sc, p)
    iz * sc.n_phi + iφ
end

# =====================================================================
# Loading
# =====================================================================

"""
    load_solid(d) -> Solid

Build a solid from a JSON dict, dispatching on its `shape` field: `"cylinder"`
(radius_cm, half_length_cm), `"box"` (half_x_cm, half_y_cm, half_z_cm),
`"cyl_shell"` (r_inner_cm, wall_thickness_cm, half_length_cm) or `"sphere"`
(radius_cm). Any other shape is rejected rather than silently mis-loaded.
"""
function load_solid(d)::Solid
    shape = get(d, "shape", "cylinder")
    if shape == "cylinder"
        Cylinder(Float64(d["radius_cm"]), Float64(d["half_length_cm"]))
    elseif shape == "box"
        Box(Float64(d["half_x_cm"]), Float64(d["half_y_cm"]), Float64(d["half_z_cm"]))
    elseif shape == "cyl_shell"
        CylShell(Float64(d["r_inner_cm"]), Float64(d["wall_thickness_cm"]),
                 Float64(d["half_length_cm"]))
    elseif shape == "sphere"
        Sphere(Float64(d["radius_cm"]))
    else
        error("unsupported solid shape '$shape' (cylinder, box, cyl_shell, sphere)")
    end
end

"""
    _load_volume(d, materials, default_name) -> PhysicalVolume

Build a placed physical volume from one geometry section `d`: the shape comes from
`load_solid`, the material name is resolved against `materials`, and the placement is
read from `position_cm` (default the origin = scanner centre).
"""
function _load_volume(d, materials::Dict{String,Material},
                      default_name::AbstractString)::PhysicalVolume
    sol = load_solid(d)
    matname = String(d["material"])
    haskey(materials, matname) ||
        error("material '$matname' not found in materials")
    mat = materials[matname]
    pos = Tuple(Float64.(get(d, "position_cm", [0.0, 0.0, 0.0])))::NTuple{3,Float64}
    PhysicalVolume(LogicalVolume(String(get(d, "name", default_name)), sol, mat), pos)
end

"Build a Scanner from its geometry section: a placed CylShell + n_phi, n_z."
function _load_scanner(d, materials::Dict{String,Material})::Scanner
    vol = _load_volume(d, materials, "scanner")
    vol.logical.solid isa CylShell ||
        error("scanner must be a 'cyl_shell' (got '$(get(d, "shape", "?"))')")
    Scanner(vol, Int(d["n_phi"]), Int(d["n_z"]))
end

"""
    Geometry

The full geometry the photons traverse. `world` is the absolute mother volume — a
cylinder of (non-interacting) Air enclosing everything; `phantom` is a daughter,
and `scanner` (the detector ring, optional) is a further daughter. This is the
volume set a multi-volume navigator will walk.
"""
# Parametric on the world/phantom solid types (the scanner volume is already a
# concrete PhysicalVolume{CylShell}) so the navigator's volume calls stay concretely
# typed — needed for the allocation-free `navigate_single_photons` hot path.
struct Geometry{W<:PhysicalVolume,P<:PhysicalVolume}
    world::W                              # the Air mother volume enclosing all daughters
    phantom::P                            # daughter
    scanner::Union{Scanner,Nothing}       # daughter (the detector ring), if present
end

"""
    load_geometry(path, materials) -> Geometry

Load the full geometry from a JSON file whose named sections each describe one
volume. The `world` (mother) and `phantom` sections are required; the `scanner`
section (the detector ring) is optional.
"""
function load_geometry(path::AbstractString,
                       materials::Dict{String,Material})::Geometry
    d = open(io -> JSON.parse(io), path, "r")
    haskey(d, "world")   || error("geometry file $path has no 'world' section")
    haskey(d, "phantom") || error("geometry file $path has no 'phantom' section")
    scanner = haskey(d, "scanner") ? _load_scanner(d["scanner"], materials) : nothing
    Geometry(_load_volume(d["world"],   materials, "world"),
             _load_volume(d["phantom"], materials, "phantom"),
             scanner)
end
