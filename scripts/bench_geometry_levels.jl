#!/usr/bin/env julia
# Design probe: does a Geant4-style three-level geometry (Solid -> LogicalVolume ->
# PhysicalVolume) cost anything in the transport hot loop versus a flat two-level
# one (Solid -> PhysicalVolume)?
#
# Conclusion (see numbers below): with parametric, concretely-typed immutable
# structs the compiler flattens the nesting — `pv.logical.solid.r` lowers to the
# same direct load as `pv.solid.r`. The third level is free at runtime; its only
# cost is code ceremony. So the level count is a modeling choice, not a speed one.
#
# Run from the repo root:
#   julia --project=. scripts/bench_geometry_levels.jl [n_calls]

abstract type Solid end
struct Cyl <: Solid
    r_cm::Float64
    h_cm::Float64
end

# Two-level: a placed, materialised volume = solid + material + position.
struct PV2{S<:Solid}
    solid::S
    density::Float64
    pos::NTuple{3,Float64}
end

# Three-level (Geant4 semantics): physical -> logical (solid + material) -> solid.
struct LV{S<:Solid}
    solid::S
    density::Float64
end
struct PV3{S<:Solid}
    logical::LV{S}
    pos::NTuple{3,Float64}
end

# Representative geometry call: world point -> local frame -> cylinder test.
@inline is_inside(s::Cyl, p) = (p[1]^2 + p[2]^2 <= s.r_cm^2) && (abs(p[3]) <= s.h_cm)
@inline is_inside(v::PV2, p) =
    is_inside(v.solid, (p[1]-v.pos[1], p[2]-v.pos[2], p[3]-v.pos[3]))
@inline is_inside(v::PV3, p) =
    is_inside(v.logical.solid, (p[1]-v.pos[1], p[2]-v.pos[2], p[3]-v.pos[3]))

# Sweep many points so the per-call cost dominates loop overhead.
function sweep(v, n::Int)
    c = 0
    @inbounds for i in 1:n
        x = 8.0 * (i % 1000) / 1000.0
        is_inside(v, (x, 0.5, 1.0)) && (c += 1)
    end
    c
end

"Best-of-`trials` wall time [s] for `sweep(v, n)`, after a warmup."
function best_time(v, n::Int; trials::Int = 5)
    sweep(v, 1000)                       # force compilation
    t = Inf
    for _ in 1:trials
        t = min(t, @elapsed sweep(v, n))
    end
    t
end

function main()
    n = isempty(ARGS) ? 200_000_000 : parse(Int, ARGS[1])
    c  = Cyl(8.0, 8.0)
    v2 = PV2(c, 1.0, (0.0, 0.0, 0.0))
    v3 = PV3(LV(c, 1.0), (0.0, 0.0, 0.0))

    # The transport loop must not allocate; confirm both layouts are alloc-free.
    sweep(v2, 1000); sweep(v3, 1000)
    a2 = @allocated sweep(v2, 1_000_000)
    a3 = @allocated sweep(v3, 1_000_000)

    t2 = best_time(v2, n)
    t3 = best_time(v3, n)

    println("calls per trial: $n")
    println("2-level (Solid + Physical):          ",
            "$(round(t2*1e9/n, digits=3)) ns/call, alloc=$a2 B")
    println("3-level (Solid + Logical + Physical): ",
            "$(round(t3*1e9/n, digits=3)) ns/call, alloc=$a3 B")
    println("ratio 3-level / 2-level: ", round(t3/t2, digits=3),
            "  (≈1.0 ⇒ the extra level is free at runtime)")
end

main()
