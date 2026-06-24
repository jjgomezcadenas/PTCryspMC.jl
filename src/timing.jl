# Timing: the per-gamma timestamp = annihilation time + time-of-flight + scintillation jitter.
# The jitter is NOT a free Gaussian σ_t: it's the arrival of the gamma's FIRST DETECTED
# scintillation photon, a physical distribution set by the crystal's light yield, decay
# mixture and photodetection efficiency (all carried on the Material). See docs/PTCryspMC_phys.tex.

"Speed of light [mm/ns] — for time-of-flight."
const C_MM_PER_NS = 299.792458

"""
    first_photon_jitter(mat, E_MeV, rng) -> Float64

Timing jitter [ns] of a gamma depositing `E_MeV` in scintillator `mat`: the arrival time of
its first detected scintillation photon. The crystal emits `N_det = light_yield·E·PDE`
detected photons, each ~Exp at the mixture's effective initial rate `r0 = Σ wₖ/τₖ`; the first
of `N_det` such draws is Exp(N_det·r0), so `jitter = −ln(u)/(N_det·r0)`, mean `1/(N_det·r0)`.
(Early-time limit — the first photon lands at ~0.1 ns ≪ τ_scint — so a single rate `r0`
handles one or two decay components uniformly.) Non-scintillators (no yield / decay / PDE)
return `0.0`. Allocation-free.
"""
@inline function first_photon_jitter(mat::Material, E_MeV::Float64, rng::AbstractRNG)::Float64
    N_det = mat.light_yield * E_MeV * mat.pde
    N_det > 0.0 || return 0.0
    r0 = 0.0                                   # Σ wₖ/τₖ over the decay mixture
    @inbounds for k in eachindex(mat.scint_decay_ns)
        τ = mat.scint_decay_ns[k]
        τ > 0.0 && (r0 += mat.scint_decay_w[k] / τ)
    end
    r0 > 0.0 || return 0.0
    # Deferred (dev/status.md #1): rand() can be exactly 0 (p≈2⁻⁵³) → -log(0) = Inf. Fix is
    # -log(1-rand), but it changes every jitter value → fold in only when singles are regenerated.
    -log(rand(rng)) / (N_det * r0)
end

"Time-of-flight [ns] from `emit` to `hit` (3-tuples, mm): ‖hit − emit‖ / c."
@inline function tof_ns(emit, hit)::Float64
    dx = Float64(hit[1]) - Float64(emit[1])
    dy = Float64(hit[2]) - Float64(emit[2])
    dz = Float64(hit[3]) - Float64(emit[3])
    sqrt(dx*dx + dy*dy + dz*dz) / C_MM_PER_NS
end

"""
    photon_timestamp(t_annih_s, emit, hit, E_MeV, mat, rng) -> Float64

The measured timestamp [ns] of one detected gamma: annihilation time + time-of-flight +
scintillation first-photon jitter. `t_annih_s` is the event annihilation time in SECONDS
(from `event_time`); it is converted to ns here so the returned stamp — and the `_ns` LOR
columns — are uniformly nanoseconds.
"""
@inline function photon_timestamp(t_annih_s::Float64, emit, hit, E_MeV::Float64,
                                  mat::Material, rng::AbstractRNG)::Float64
    t_annih_s * 1.0e9 + tof_ns(emit, hit) + first_photon_jitter(mat, E_MeV, rng)
end
