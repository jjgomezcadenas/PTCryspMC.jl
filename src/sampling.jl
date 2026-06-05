# Monte Carlo samplers. Adapted from LXeMC (src/sampling.jl), photon-only.

const ME = 0.51099895  # electron rest mass [MeV]

"Free-flight distance [cm] from an exponential with macroscopic cross section Σ [cm^-1]."
sample_distance(Σ::Float64, rng::AbstractRNG)::Float64 = -log(rand(rng)) / Σ

"Pick :compton / :pair / :photoelectric from the three branching probabilities."
function sample_process(P_C::Float64, P_P::Float64, P_Ph::Float64, rng::AbstractRNG)::Symbol
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
    abs(rd[3]) > 0.99999 && return local_vec .* sign(rd[3])
    rp = sqrt(rd[1]^2 + rd[2]^2)
    e1 = Float64[-rd[3]*rd[1], -rd[3]*rd[2], -(rd[3]^2 - 1.0)] ./ rp
    e2 = Float64[-rd[2], rd[1], 0.0] ./ rp
    local_vec[1] .* e1 .+ local_vec[2] .* e2 .+ local_vec[3] .* rd
end
