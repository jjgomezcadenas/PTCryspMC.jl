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
    t_far > 0.0 ? t_far : Inf
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
# Loading
# =====================================================================

"""
    load_solid(d) -> Solid

Build a solid from a JSON dict, dispatching on its `shape` field. Only
`"cylinder"` (radius_cm, half_length_cm) is supported so far — any other shape
is rejected rather than silently mis-loaded.
"""
function load_solid(d)::Solid
    shape = get(d, "shape", "cylinder")
    if shape == "cylinder"
        Cylinder(Float64(d["radius_cm"]), Float64(d["half_length_cm"]))
    elseif shape == "box"
        Box(Float64(d["half_x_cm"]), Float64(d["half_y_cm"]), Float64(d["half_z_cm"]))
    else
        error("unsupported solid shape '$shape' (cylinder, box)")
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

"""
    Geometry

The full geometry — the world the photons traverse. For now it holds only the
phantom; the detector ring is added as a second component later. This is the seed
of the volume list a multi-volume navigator will walk.
"""
struct Geometry
    phantom::PhysicalVolume
end

"""
    load_geometry(path, materials) -> Geometry

Load the full geometry from a JSON file whose named sections each describe one
component. So far only the `phantom` section is read; the detector ring will be a
second section.
"""
function load_geometry(path::AbstractString,
                       materials::Dict{String,Material})::Geometry
    d = open(io -> JSON.parse(io), path, "r")
    haskey(d, "phantom") || error("geometry file $path has no 'phantom' section")
    Geometry(_load_volume(d["phantom"], materials, "phantom"))
end
