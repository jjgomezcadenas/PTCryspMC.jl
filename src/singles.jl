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

# =====================================================================
# Integer quantization (the compact singles encoding)
# =====================================================================
# Positions are stored in 0.1 mm units and energy in 0.1 keV units as integers — far finer
# than the ~2 mm / ~50 keV detector resolution (and downstream smearing dwarfs the grid), so
# the quantization is lossless in practice. Int16 @ 0.1 spans ±3276.7 mm / 0–3276.7 keV, with
# a bounds check so an oversized geometry fails loud instead of silently wrapping. The scale
# is stored in the file (HDF5 attribute) so readers decode by multiplying.

const XYZ_SCALE_MM = 0.1     # mm per integer unit  → int = round(mm / 0.1)
const E_SCALE_KEV  = 0.1     # keV per integer unit
const _I16MAX = 32767

@inline function encode_xyz_mm(mm::Real)::Int16
    u = round(Int, mm / XYZ_SCALE_MM)
    abs(u) <= _I16MAX || error("position $mm mm exceeds Int16 @ $XYZ_SCALE_MM mm " *
                               "(±$(_I16MAX*XYZ_SCALE_MM) mm) — widen the scale")
    Int16(u)
end
@inline function encode_e_keV(keV::Real)::Int16
    u = round(Int, keV / E_SCALE_KEV)
    (0 <= u <= _I16MAX) || error("energy $keV keV exceeds Int16 @ $E_SCALE_KEV keV")
    Int16(u)
end
@inline decode_xyz(i::Integer)::Float64 = i * XYZ_SCALE_MM   # → mm
@inline decode_e(i::Integer)::Float64   = i * E_SCALE_KEV    # → keV

# =====================================================================
# SinglesBuffer: the columnar, quantized in-memory store + binary part I/O
# =====================================================================
# One thread's chunk accumulates here, then dumps a binary part-file; the serial HDF5 packer
# reads parts back. Plain file I/O (no HDF5 in the parallel region — the HDF5 C library is not
# thread-safe). Columns are written/read in the canonical order of `singles_columns`.

struct SinglesBuffer
    event::Vector{Int32}
    gamma::Vector{Int8}
    x::Vector{Int16}; y::Vector{Int16}; z::Vector{Int16}
    e::Vector{Int16}
    iz::Vector{Int16}; iphi::Vector{Int16}
    nblocks::Vector{Int8}
    phantom_scatter::Vector{Int8}
    x0::Vector{Int16}; y0::Vector{Int16}; z0::Vector{Int16}
end
SinglesBuffer() = SinglesBuffer(Int32[], Int8[], Int16[], Int16[], Int16[], Int16[],
                                Int16[], Int16[], Int8[], Int8[], Int16[], Int16[], Int16[])
Base.length(b::SinglesBuffer) = length(b.event)

"The 13 singles columns in canonical order: (name, vector). The schema + I/O ordering."
singles_columns(b::SinglesBuffer) = (
    ("event_number", b.event), ("gamma", b.gamma),
    ("x_mm", b.x), ("y_mm", b.y), ("z_mm", b.z), ("e_keV", b.e),
    ("iz", b.iz), ("iphi", b.iphi), ("nblocks", b.nblocks), ("phantom_scatter", b.phantom_scatter),
    ("x0_mm", b.x0), ("y0_mm", b.y0), ("z0_mm", b.z0))

"Append one detected photon (the `navigate_single_photons` summary `s` + emission point `pos0`, both cm), quantized."
function push_single!(b::SinglesBuffer, ev::Integer, g::Integer, s, pos0)
    push!(b.event, Int32(ev)); push!(b.gamma, Int8(g))
    push!(b.x, encode_xyz_mm(s.x * 10)); push!(b.y, encode_xyz_mm(s.y * 10)); push!(b.z, encode_xyz_mm(s.z * 10))
    push!(b.e, encode_e_keV(s.e * 1000))
    push!(b.iz, Int16(s.iz)); push!(b.iphi, Int16(s.iphi))
    push!(b.nblocks, Int8(s.nblocks)); push!(b.phantom_scatter, Int8(s.phscat ? 1 : 0))
    push!(b.x0, encode_xyz_mm(pos0[1] * 10)); push!(b.y0, encode_xyz_mm(pos0[2] * 10)); push!(b.z0, encode_xyz_mm(pos0[3] * 10))
    b
end

"Dump the buffer's columns as raw bytes (canonical order) to a binary part-file."
write_part(io::IO, b::SinglesBuffer) = foreach(c -> write(io, c[2]), singles_columns(b))

"Read `n` rows from a binary part-file (canonical column order) into a fresh `SinglesBuffer`."
function read_part(io::IO, n::Int)::SinglesBuffer
    rd(T) = read!(io, Vector{T}(undef, n))
    SinglesBuffer(rd(Int32), rd(Int8), rd(Int16), rd(Int16), rd(Int16), rd(Int16),
                  rd(Int16), rd(Int16), rd(Int8), rd(Int8), rd(Int16), rd(Int16), rd(Int16))
end
