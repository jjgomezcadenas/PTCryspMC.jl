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
    sigma_xyz::Float64       # mm (position model 1: single Gaussian σ per axis)
    eres::Float64            # fractional FWHM at 511 keV
    emin::Float64            # minimum gamma energy [keV] (0 = off)
    apply_window::Bool
    win_half::Float64        # symmetric window half-width about 511 keV [keV]
    pos_model::Int           # 1 = single Gaussian (σ_xyz), 2 = core/tail mixture (pos2)
    mix_s1::Float64          # model 2: core Gaussian σ [mm]
    mix_s2::Float64          # model 2: tail Gaussian σ [mm]
    mix_f::Float64           # model 2: core fraction (tier-resolved)
    sel_eff::Float64         # per-gamma selection efficiency (1 = keep all; the model-2 tier)
end

# Back-compat: the 5-arg form is position model 1, no selection (all existing call sites).
Response(sigma_xyz, eres, emin, apply_window, win_half) =
    Response(sigma_xyz, eres, emin, apply_window, win_half, 1, 0.0, 0.0, 0.0, 1.0)

"Smear a hit position by the response's position model (1: single Gaussian; 2: core/tail mixture)."
@inline smear_hit(p, r::Response, rng) =
    r.pos_model == 2 ? smear_position_mix(p, r.mix_s1, r.mix_s2, r.mix_f, rng) :
                       smear_position(p, r.sigma_xyz, rng)

# A (smeared) gamma energy passes the energy selection: above the minimum AND, if a symmetric
# window is set, within it.
@inline pass_energy(e::Float64, r::Response) =
    (r.emin <= 0.0 || e >= r.emin) && (!r.apply_window || abs(e - 511.0) <= r.win_half)

# Per-gamma accumulator, reset at each event boundary. Holds ONLY what the fill_* calls
# accumulate: the first scanner deposit fixes the LOR point and block; later deposits add energy
# and, in a different block, flag overspill. The two non-accumulable per-gamma scalars — the
# phantom-scatter count and the timestamp `t` (which needs the gamma's *total* energy) — are NOT
# stored here; the caller passes them to `finish_event!`. So there is no half-filled state: a
# `GammaAcc` is valid the moment filling stops.
mutable struct GammaAcc
    reached::Bool
    x::Float64; y::Float64; z::Float64
    iz::Int; iphi::Int
    e::Float64
    overspill::Bool          # a deposit landed in a second, different block
end
GammaAcc() = GammaAcc(false, 0.0, 0.0, 0.0, -1, -1, 0.0, false)

function reset!(a::GammaAcc)
    a.reached = false; a.x = a.y = a.z = 0.0; a.iz = a.iphi = -1
    a.e = 0.0; a.overspill = false
    a
end

"A gamma passes if it reached the detector and stayed within one crystal block."
contained_one(a::GammaAcc) = a.reached && !a.overspill

"""
    fill_full!(a, x, y, z, e_dep, iz, iphi)

Accumulate one FULL-stack scanner deposit into `a`: the first deposit fixes the LOR point and
block; a later deposit in a different block flags overspill; energy sums. The phantom-scatter
count and the timestamp are NOT set here — the caller tracks them per gamma (nscat summed over the
phantom rows; `t` stamped once the energy is complete) and passes them to `finish_event!`.
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
    fill_singles!(a, x, y, z, e, iz, iphi, nblocks)

Fill `a` directly from a SINGLES row — the row already is the formed hit (`reached`, LOR point,
summed energy, `overspill = nblocks > 1`). The per-gamma phantom-scatter count and timestamp
`t_rel` are NOT stored on the accumulator; the caller carries them and passes them to
`finish_event!`.
"""
@inline function fill_singles!(a::GammaAcc, x::Float64, y::Float64, z::Float64, e::Float64,
                               iz::Int, iphi::Int, nblocks::Int)
    a.reached = true; a.x = x; a.y = y; a.z = z; a.iz = iz; a.iphi = iphi; a.e = e
    a.overspill = nblocks > 1
    a
end

# Truth codes carried in the LOR list.
const TRUTH_TRUE    = Int8(0)
const TRUTH_SCATTER = Int8(1)
const TRUTH_RANDOM  = Int8(2)     # set by the randoms pass (build_randoms → reco_lors)

"""
    finish_event!(emit, ev, g1, g2, t1, t2, nscat1, nscat2, x0, y0, z0, resp, rng) -> (emitted, is_true)

Finalise one annihilation's two gammas into a LOR: require both contained in one block, smear
energy + position, apply the energy selection, and on acceptance call
`emit(ev, x1,y1,z1,e1,t1,iz1,iphi1,nscat1, x2,y2,z2,e2,t2,iz2,iphi2,nscat2, dt, x0,y0,z0, truth)`.

The per-gamma timestamps `t1,t2` (relative to the decay = `TOF + jitter`) and phantom-scatter
counts `nscat1,nscat2` are passed in, NOT carried on the accumulator — `GammaAcc` holds only what
the `fill_*` calls accumulate (geometry + energy + overspill), so there is no half-filled state to
forget. `truth = TRUTH_TRUE` unless either gamma scattered in the phantom (`nscatᵢ > 0`); the counts
also separate single (`nscat1+nscat2 == 1`) from multiple (`≥ 2`) scatter. Smearing does not change
the tag. `dt = (t1−t2) − TOF_diff` is the TOF-corrected timing residual (the signed jitter
difference, centred at 0; `TOF_diff` = the geometric time-of-flight difference from the common
emission point). Returns whether a LOR was emitted and whether it is a true coincidence.
"""
function finish_event!(emit, ev::Int, g1::GammaAcc, g2::GammaAcc,
                       t1::Float64, t2::Float64, nscat1::Int, nscat2::Int,
                       x0::Float64, y0::Float64, z0::Float64, resp::Response, rng)
    (contained_one(g1) && contained_one(g2)) || return (false, false)
    if resp.sel_eff < 1.0        # model-2 selection tier: per-gamma Bernoulli (guarded so
        k1 = rand(rng) < resp.sel_eff   # sel_eff = 1 draws nothing — legacy streams unchanged)
        k2 = rand(rng) < resp.sel_eff
        (k1 && k2) || return (false, false)
    end
    e1 = smear_energy(g1.e, resp.eres, rng)
    e2 = smear_energy(g2.e, resp.eres, rng)
    (pass_energy(e1, resp) && pass_energy(e2, resp)) || return (false, false)
    x1, y1, z1 = smear_hit((g1.x, g1.y, g1.z), resp, rng)
    x2, y2, z2 = smear_hit((g2.x, g2.y, g2.z), resp, rng)
    is_true = (nscat1 == 0 && nscat2 == 0)         # true iff neither gamma scattered in the phantom
    truth = is_true ? TRUTH_TRUE : TRUTH_SCATTER   # nscat1+nscat2: 1 = single, ≥2 = multiple scatter

    # The residual subtracts the (signed) geometric TOF difference from the (signed) timestamp
    # difference, so the geometry cancels exactly and dt = (t1−t2) − TOF_diff is the pure jitter
    # difference (signed, centred at 0). The abs must NOT wrap the timestamp difference — that would
    # leave ~2·TOF_diff of geometry for half the pairs.
    emit_pt = (x0, y0, z0)
    tof_diff = tof_ns(emit_pt, (g1.x, g1.y, g1.z)) - tof_ns(emit_pt, (g2.x, g2.y, g2.z))
    dt = (t1 - t2) - tof_diff

    emit(ev, x1, y1, z1, e1, t1, g1.iz, g1.iphi, nscat1,
             x2, y2, z2, e2, t2, g2.iz, g2.iphi, nscat2, dt, x0, y0, z0, truth)
    (true, is_true)
end
