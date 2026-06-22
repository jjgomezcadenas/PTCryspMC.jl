# Coincidence (LOR) selection: pair the two photons of an annihilation into a line of response.
# The shared core behind both LOR builders — build_coincidences.jl (full stack → CSV) and
# build_true_coincidences_from_singles.jl (singles → LOR HDF5). Sink-agnostic: `finish_event!` emits
# each accepted LOR via a callback (a CSV/HDF5 row in the drivers, a vector push in the tests),
# so the selection physics lives in one place.
#
# A coincidence is accepted when both gammas reached the detector and each stayed within ONE
# crystal block (no overspill); partial energy is allowed (an escaping gamma from a single
# crystal still gives a clean LOR point). In detector mode the hits are smeared (σ_xyz, FWHM(E))
# and an energy selection (symmetric window or `emin` minimum) is applied to the smeared energy.

# Detector response. All-off = truth mode (no smearing, no energy cut).
struct Response
    sigma_xyz::Float64       # mm
    eres::Float64            # fractional FWHM at 511 keV
    emin::Float64            # minimum gamma energy [keV] (0 = off)
    apply_window::Bool
    win_half::Float64        # symmetric window half-width about 511 keV [keV]
end

# A (smeared) gamma energy passes the energy selection: above the minimum AND, if a symmetric
# window is set, within it.
@inline pass_energy(e::Float64, r::Response) =
    (r.emin <= 0.0 || e >= r.emin) && (!r.apply_window || abs(e - 511.0) <= r.win_half)

# Per-gamma accumulator, reset at each event boundary. The first scanner deposit fixes the LOR
# point and block; later deposits add energy and, in a different block, flag overspill.
mutable struct GammaAcc
    reached::Bool
    x::Float64; y::Float64; z::Float64
    iz::Int; iphi::Int
    e::Float64
    overspill::Bool          # a deposit landed in a second, different block
    phscat::Bool
end
GammaAcc() = GammaAcc(false, 0.0, 0.0, 0.0, -1, -1, 0.0, false, false)

function reset!(a::GammaAcc)
    a.reached = false; a.x = a.y = a.z = 0.0; a.iz = a.iphi = -1
    a.e = 0.0; a.overspill = false; a.phscat = false
    a
end

"A gamma passes if it reached the detector and stayed within one crystal block."
contained_one(a::GammaAcc) = a.reached && !a.overspill

"""
    fill_full!(a, x, y, z, e_dep, iz, iphi)

Accumulate one FULL-stack scanner deposit into `a`: the first deposit fixes the LOR point and
block; a later deposit in a different block flags overspill; energy sums. (The phantom-scatter
flag is per-photon and set by the caller from the row, not here.)
"""
@inline function fill_full!(a::GammaAcc, x::Float64, y::Float64, z::Float64,
                            e_dep::Float64, iz::Int, iphi::Int)
    if !a.reached
        a.reached = true; a.x = x; a.y = y; a.z = z; a.iz = iz; a.iphi = iphi
    elseif iz != a.iz || iphi != a.iphi
        a.overspill = true
    end
    a.e += e_dep
    a
end

"""
    fill_singles!(a, x, y, z, e, iz, iphi, nblocks, phscat)

Fill `a` directly from a SINGLES row — the row already is the formed hit (`reached`, LOR point,
summed energy, `overspill = nblocks > 1`, phantom-scatter flag).
"""
@inline function fill_singles!(a::GammaAcc, x::Float64, y::Float64, z::Float64, e::Float64,
                               iz::Int, iphi::Int, nblocks::Int, phscat::Bool)
    a.reached = true; a.x = x; a.y = y; a.z = z; a.iz = iz; a.iphi = iphi; a.e = e
    a.overspill = nblocks > 1; a.phscat = phscat
    a
end

# Truth codes carried in the LOR list.
const TRUTH_TRUE    = Int8(0)
const TRUTH_SCATTER = Int8(1)
const TRUTH_RANDOM  = Int8(2)     # reserved for the Step-5 randoms pass

"""
    finish_event!(emit, ev, g1, g2, x0, y0, z0, resp, rng, timing) -> (emitted, is_true)

Finalise one annihilation's two gammas into a LOR: require both contained in one block, smear
energy + position, apply the energy selection, time-stamp each gamma, and on acceptance call
`emit(ev, x1,y1,z1,e1,t1,iz1,iphi1, x2,y2,z2,e2,t2,iz2,iphi2, dt, x0,y0,z0, truth)`,
`truth` ∈ `TRUTH_TRUE`/`TRUTH_SCATTER` (smearing does not change the tag).

`timing::EventTiming` supplies the crystal (its scintillation jitter). Timestamps `t1,t2` are
**relative to the decay** — `TOF + jitter`, from the **true** hit positions and **true** deposited
energies (the shared annihilation time cancels in a same-pair difference, so it is dropped;
absolute = `event_time(timing.act, ev) + t`). The residual `dt = |t1−t2| − TOF_diff` is the
timing-resolution residual (`TOF_diff` = geometric time-of-flight difference). Returns whether a
LOR was emitted and whether it is a true coincidence.
"""
function finish_event!(emit, ev::Int, g1::GammaAcc, g2::GammaAcc,
                       x0::Float64, y0::Float64, z0::Float64, resp::Response, rng,
                       timing)   # an EventTiming (untyped: this file is included before timing.jl)
    (contained_one(g1) && contained_one(g2)) || return (false, false)
    e1 = smear_energy(g1.e, resp.eres, rng)
    e2 = smear_energy(g2.e, resp.eres, rng)
    (pass_energy(e1, resp) && pass_energy(e2, resp)) || return (false, false)
    x1, y1, z1 = smear_position((g1.x, g1.y, g1.z), resp.sigma_xyz, rng)
    x2, y2, z2 = smear_position((g2.x, g2.y, g2.z), resp.sigma_xyz, rng)
    is_true = !(g1.phscat || g2.phscat)
    truth = is_true ? TRUTH_TRUE : TRUTH_SCATTER

    # Per-gamma timestamp RELATIVE to the decay: TOF (true hit ← emit) + scintillation jitter (from
    # true deposited energy, keV→MeV). The annihilation time t_annih is common to both photons of a
    # pair, so it cancels in any same-pair difference and is dropped here — keeping t1,t2 at the ns
    # scale (Float32-exact). Absolute time = event_time(timing.act, ev) + t_rel, recomputable from
    # the stored event_number. DT subtracts the geometric TOF difference → the timing-res residual.
    emit_pt = (x0, y0, z0)
    tof1 = tof_ns(emit_pt, (g1.x, g1.y, g1.z))
    tof2 = tof_ns(emit_pt, (g2.x, g2.y, g2.z))
    t1 = tof1 + first_photon_jitter(timing.mat, g1.e * 1e-3, rng)
    t2 = tof2 + first_photon_jitter(timing.mat, g2.e * 1e-3, rng)
    dt = abs(t1 - t2) - (tof1 - tof2)

    emit(ev, x1, y1, z1, e1, t1, g1.iz, g1.iphi,
             x2, y2, z2, e2, t2, g2.iz, g2.iphi, dt, x0, y0, z0, truth)
    (true, is_true)
end
