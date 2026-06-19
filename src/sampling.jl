# Monte Carlo samplers, photon-only: free path, process choice, Compton scatter.

const ME = 0.51099895  # electron rest mass [MeV]

"Free-flight distance [cm] from an exponential with macroscopic cross section Σ [cm^-1]."
sample_distance(Σ::Float64, rng::AbstractRNG)::Float64 = -log(rand(rng)) / Σ

"""
Pick :compton / :photoelectric / :pair from the branching probabilities, in the
(C, Ph, P) order returned by `sigma_macro`. Photoelectric is the catch-all so
that a zero-width pair bucket (P below threshold) can never absorb a rounding
leftover in `r`.
"""
function sample_process(P_C::Float64, P_Ph::Float64, P_P::Float64, rng::AbstractRNG)::Symbol
    r = rand(rng)
    r < P_C && return :compton
    r < P_C + P_P && return :pair
    :photoelectric
end

"""
    sample_compton(E, rng) -> (E_scattered, cos_theta)

Klein-Nishina Compton scattering of a photon of energy `E` [MeV], by the
Butcher-Messel composition+rejection method (as in Geant4).
"""
function sample_compton(E::Float64, rng::AbstractRNG)::Tuple{Float64,Float64}
    eps0 = ME / (ME + 2.0 * E)
    a1 = -log(eps0)
    a2 = (1.0 - eps0^2) / 2.0
    while true
        r1 = rand(rng); r2 = rand(rng); r3 = rand(rng)
        eps = r1 < a1 / (a1 + a2) ? eps0^r2 : sqrt(eps0^2 + (1.0 - eps0^2) * r2)
        one_minus_cos = (1.0 - eps) * ME / (E * eps)
        sin2 = clamp(one_minus_cos * (2.0 - one_minus_cos), 0.0, 1.0)
        g = 1.0 - eps * sin2 / (1.0 + eps^2)
        g >= r3 && return (eps * E, 1.0 - one_minus_cos)
    end
end

"""
    rotate_to_global(local_vec, ref_dir) -> Vector{Float64}

Rotate a unit vector from a local frame (z = ref_dir) into the global frame.
"""
function rotate_to_global(local_vec, ref_dir)::Vector{Float64}
    n = sqrt(ref_dir[1]^2 + ref_dir[2]^2 + ref_dir[3]^2)
    rd = ref_dir ./ n
    # Near-pole: the local frame ≈ the global frame (up to a flip). Build a Vector
    # explicitly so a tuple `local_vec` (e.g. from the acollinearity tilt) still returns
    # the declared Vector{Float64}, not a tuple.
    if abs(rd[3]) > 0.99999
        s = sign(rd[3])
        return Float64[local_vec[1]*s, local_vec[2]*s, local_vec[3]*s]
    end
    rp = sqrt(rd[1]^2 + rd[2]^2)
    e1 = Float64[-rd[3]*rd[1], -rd[3]*rd[2], -(rd[3]^2 - 1.0)] ./ rp
    e2 = Float64[-rd[2], rd[1], 0.0] ./ rp
    local_vec[1] .* e1 .+ local_vec[2] .* e2 .+ local_vec[3] .* rd
end

"""
    sample_interaction(E, dir, ΣC, ΣPh, ΣP, rng) -> (process, e_dep, new_dir, new_E)

The physics at a single interaction point, record-type agnostic. Given the photon
energy `E` [MeV], its direction, and the already-computed macroscopic cross sections
(so `sigma_macro` is called once per step in the caller's loop), pick the process and
return what is deposited and the photon's new state:

- `:compton`  — deposit the electron recoil `E − E'`, photon continues as `(new_dir, E')`;
- `:photoelectric` / `:pair` — full absorption, `new_E = 0` (`new_dir` returned unchanged).

Shared by the single-volume `propagate_photon` and the multi-volume `navigate_photon`;
the below-cut / stop bookkeeping is a loop concern and stays in each caller.
"""
function sample_interaction(E::Float64, dir, ΣC::Float64, ΣPh::Float64, ΣP::Float64,
                            rng::AbstractRNG)
    Σ = ΣC + ΣPh + ΣP
    proc = sample_process(ΣC / Σ, ΣPh / Σ, ΣP / Σ, rng)
    if proc === :compton
        Eprime, cosθ = sample_compton(E, rng)
        ϕ = 2π * rand(rng)
        sinθ = sqrt(max(0.0, 1.0 - cosθ^2))
        ndir = rotate_to_global(Float64[sinθ*cos(ϕ), sinθ*sin(ϕ), cosθ], dir)
        return (:compton, E - Eprime, ndir, Eprime)
    else
        return (proc, E, dir, 0.0)        # photoelectric / pair: full absorption
    end
end
