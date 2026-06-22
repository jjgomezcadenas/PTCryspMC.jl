# Randoms: accidental (cross-decay) coincidences. Two singles from DIFFERENT annihilations whose
# ABSOLUTE arrival times fall within the coincidence window τ form a false LOR. Absolute time
# (= event_time(decay) + t_rel) is what decides overlap — decays are spread over the whole
# acquisition (~600 s) while τ is ~ns, so only singles whose decays nearly coincide can pair.
# This is the one place the absolute clock is needed; it stays in memory — nothing absolute is
# stored (the LOR keeps small decay-relative times).

"""
    pair_randoms(emit, t_abs, ev, tau) -> npairs

Sort singles by absolute time `t_abs` [ns] and emit each CROSS-event pair within the coincidence
window `tau` [ns]: `emit(i, j, Δ)` with `i` the earlier single, `j` the later, and `Δ = t_abs[j] −
t_abs[i] ≤ tau` (the pair's `|Δt|`, also the LOR's `dT`). Same-event pairs (the two gammas of one
decay) are skipped — those are the trues, not randoms. Every cross-event pair in the window is
counted (multiples included). An O(n log n) sort + a forward sliding window that breaks as soon as
`Δ > tau` (the window holds ~S·τ singles, ≪ n at realistic rates). `i,j` index the ORIGINAL order.
"""
function pair_randoms(emit, t_abs::AbstractVector{Float64}, ev::AbstractVector{<:Integer},
                      tau::Float64)::Int
    order = sortperm(t_abs)
    n = length(order)
    np = 0
    @inbounds for a in 1:n
        i = order[a]; ti = t_abs[i]
        for b in (a + 1):n
            j = order[b]
            Δ = t_abs[j] - ti
            Δ > tau && break                 # sorted by time → no later single can be within τ
            ev[i] == ev[j] && continue        # same decay → a true/scatter pair, not a random
            emit(i, j, Δ)
            np += 1
        end
    end
    np
end
