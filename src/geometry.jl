# Geometry: a cylinder solid + ray-cylinder distance, and the phantom loader.
# Adapted from LXeMC (src/geometry_core.jl): the Cyl primitive and its
# distance_to_exit / distance_to_entry, here for a single positioned cylinder.

"A cylinder, axis along z, centred at `position` [cm]."
struct Cylinder
    radius_cm::Float64
    half_height_cm::Float64
    position::NTuple{3,Float64}
end
Cylinder(r::Real, h::Real) = Cylinder(Float64(r), Float64(h), (0.0, 0.0, 0.0))

function is_inside(c::Cylinder, p)::Bool
    dx = p[1] - c.position[1]; dy = p[2] - c.position[2]; dz = p[3] - c.position[3]
    (dx^2 + dy^2 <= c.radius_cm^2) && (abs(dz) <= c.half_height_cm)
end

"Distance from `pos` along `dir` to where the ray leaves the cylinder [cm] (Inf if it never does)."
function distance_to_exit(pos, dir, c::Cylinder)::Float64
    R = c.radius_cm; H = c.half_height_cm
    dx = pos[1] - c.position[1]; dy = pos[2] - c.position[2]; dz = pos[3] - c.position[3]
    t_min = Inf
    a = dir[1]^2 + dir[2]^2
    if a > 1e-20
        b = 2.0 * (dx * dir[1] + dy * dir[2])
        cc = dx^2 + dy^2 - R^2
        disc = b^2 - 4.0 * a * cc
        if disc >= 0.0
            sq = sqrt(disc)
            for t in ((-b + sq) / (2a), (-b - sq) / (2a))
                (t > 1e-10 && abs(dz + t * dir[3]) < H) && (t_min = min(t_min, t))
            end
        end
    end
    if abs(dir[3]) > 1e-20
        for zf in (H, -H)
            t = (zf - dz) / dir[3]
            if t > 1e-10
                rx = dx + t * dir[1]; ry = dy + t * dir[2]
                (rx^2 + ry^2 < R^2) && (t_min = min(t_min, t))
            end
        end
    end
    t_min
end

"A phantom = a cylinder + the name of its material."
struct Phantom
    cyl::Cylinder
    material::String
end

"Load a phantom from a JSON file (shape=cylinder; radius_cm, half_length_cm, material)."
function load_phantom(path::AbstractString)::Phantom
    d = open(path, "r") do io
        JSON.parse(io)
    end
    Phantom(Cylinder(Float64(d["radius_cm"]), Float64(d["half_length_cm"])),
            String(d["material"]))
end
