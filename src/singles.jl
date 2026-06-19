# Singles generation core: the chunked, sink-agnostic unit of work behind the production
# multi-threaded driver (scripts/simulate_source_mt.jl). Kept in the package (not the script)
# so the driver and the tests run the SAME code — the chunk partition, per-chunk RNG use and
# transport→reduce are tested directly, not re-inspected.

"""
    chunk_ranges(nevents, nchunks) -> Vector{UnitRange{Int}}

Partition events `1:nevents` into `nchunks` CONTIGUOUS ranges that exactly tile `1:nevents`
(no gap, no overlap), balanced by integer division. `nchunks` is clamped to `[1, nevents]`
(never more chunks than events). Contiguity is what keeps the output event-ordered once the
per-chunk part-files are concatenated in chunk order.
"""
function chunk_ranges(nevents::Int, nchunks::Int)::Vector{UnitRange{Int}}
    nevents >= 1 || error("nevents must be ≥ 1")
    k = max(1, min(nchunks, nevents))
    [(div((c - 1) * nevents, k) + 1):(div(c * nevents, k)) for c in 1:k]
end

"""
    singles_chunk!(emit, geom, src, E0, cut_MeV, acol, range, rng) -> nrows

Generate the singles of one event `range`: for each event draw a back-to-back pair from `src`
(acollinearity `acol`° FWHM) and navigate both photons with `navigate_single_photons` (the
allocation-free reducer), calling `emit(ev, gamma, summary, pos0)` once per photon that
reached the ring — `summary` is the `navigate_single_photons` NamedTuple, `pos0` the emission
point [cm]. The `emit` sink (a file writer in the driver, a vector push in the tests) decides
where rows go; the physics here is identical either way. Pure given `(geom, src, params, rng)`
— no threads, no I/O — so it is the unit a thread runs and the unit the tests check. Returns
the number of rows emitted.
"""
function singles_chunk!(emit, geom::Geometry, src::Source, E0::Float64, cut_MeV::Float64,
                        acol::Float64, range::UnitRange{Int}, rng::AbstractRNG)::Int
    nr = 0
    for ev in range
        pos0, d1, d2 = emit_pair(src, rng; acol_fwhm_deg=acol)
        for (g, dir) in ((1, d1), (2, d2))
            s = navigate_single_photons(geom, E0, pos0, dir, rng; egamma_cut=cut_MeV)
            s.reached || continue
            emit(ev, g, s, pos0)
            nr += 1
        end
    end
    nr
end
