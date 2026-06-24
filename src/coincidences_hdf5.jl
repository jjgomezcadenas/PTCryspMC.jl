# HDF5 container for the LOR (coincidence) list: a quantized columnar buffer + a streaming
# writer over extensible datasets (the total is unknown until the selection stream ends), and a
# slice reader. Same quantization as the singles (0.1 mm / 0.1 keV Int16, the scales stored as
# attributes), shuffle+deflate-4. The LOR values are SMEARED (detector response), so a Gaussian
# tail can fall out of range — encode with CLAMP (not the bounds-error used for singles truth).
#
# Timestamps t1_ns/t2_ns are Float32 and stored RELATIVE TO THE DECAY (= TOF + scintillation
# jitter, ~ns) — the common annihilation time is dropped so ns differences survive Float32.
# Absolute time = event_time(activity, event) + t; the activity attrs (t0_s/half_life_s/time_seed)
# are written alongside. dt_ns is the residual (t1−t2) − TOF_diff (signed, the timing resolution).

using HDF5

const COINC_BLOCK = 1 << 20     # rows buffered before a block is appended to the datasets

@inline _enc_xyz_c(mm::Real)::Int16 = Int16(clamp(round(Int, mm / XYZ_SCALE_MM), -32767, 32767))
@inline _enc_e_c(keV::Real)::Int16  = Int16(clamp(round(Int, keV / E_SCALE_KEV), 0, 32767))

struct CoincidenceBuffer
    event::Vector{Int32}; truth::Vector{Int8}
    x1::Vector{Int16}; y1::Vector{Int16}; z1::Vector{Int16}; e1::Vector{Int16}
    t1::Vector{Float32}; iz1::Vector{Int16}; iphi1::Vector{Int16}
    x2::Vector{Int16}; y2::Vector{Int16}; z2::Vector{Int16}; e2::Vector{Int16}
    t2::Vector{Float32}; iz2::Vector{Int16}; iphi2::Vector{Int16}
    dt::Vector{Float32}                       # per-pair timing residual (t1−t2) − TOF_diff [ns], signed
    x0::Vector{Int16}; y0::Vector{Int16}; z0::Vector{Int16}
end
CoincidenceBuffer() = CoincidenceBuffer(Int32[], Int8[], Int16[], Int16[], Int16[], Int16[],
    Float32[], Int16[], Int16[], Int16[], Int16[], Int16[], Int16[], Float32[], Int16[], Int16[],
    Float32[], Int16[], Int16[], Int16[])
Base.length(b::CoincidenceBuffer) = length(b.event)

"The LOR columns in canonical order (name, vector) — the schema + dataset/IO ordering."
coinc_columns(b::CoincidenceBuffer) = (
    ("event", b.event), ("truth", b.truth),
    ("x1_mm", b.x1), ("y1_mm", b.y1), ("z1_mm", b.z1), ("e1_keV", b.e1), ("t1_ns", b.t1), ("iz1", b.iz1), ("iphi1", b.iphi1),
    ("x2_mm", b.x2), ("y2_mm", b.y2), ("z2_mm", b.z2), ("e2_keV", b.e2), ("t2_ns", b.t2), ("iz2", b.iz2), ("iphi2", b.iphi2),
    ("dt_ns", b.dt),
    ("x0_mm", b.x0), ("y0_mm", b.y0), ("z0_mm", b.z0))

Base.empty!(b::CoincidenceBuffer) = (foreach(c -> empty!(c[2]), coinc_columns(b)); b)

"Append one accepted LOR (the `finish_event!` emit args), quantized."
function push_coincidence!(b::CoincidenceBuffer, ev, x1, y1, z1, e1, t1, iz1, iphi1,
                           x2, y2, z2, e2, t2, iz2, iphi2, dt, x0, y0, z0, truth)
    push!(b.event, Int32(ev)); push!(b.truth, Int8(truth))
    push!(b.x1, _enc_xyz_c(x1)); push!(b.y1, _enc_xyz_c(y1)); push!(b.z1, _enc_xyz_c(z1))
    push!(b.e1, _enc_e_c(e1)); push!(b.t1, Float32(t1)); push!(b.iz1, Int16(iz1)); push!(b.iphi1, Int16(iphi1))
    push!(b.x2, _enc_xyz_c(x2)); push!(b.y2, _enc_xyz_c(y2)); push!(b.z2, _enc_xyz_c(z2))
    push!(b.e2, _enc_e_c(e2)); push!(b.t2, Float32(t2)); push!(b.iz2, Int16(iz2)); push!(b.iphi2, Int16(iphi2))
    push!(b.dt, Float32(dt))
    push!(b.x0, _enc_xyz_c(x0)); push!(b.y0, _enc_xyz_c(y0)); push!(b.z0, _enc_xyz_c(z0))
    b
end

# Streaming HDF5 writer: buffer rows, append a block to the (extensible) datasets when full.
mutable struct CoincidenceWriter
    file::HDF5.File
    dsets::Dict{String,HDF5.Dataset}
    buf::CoincidenceBuffer
    n::Int
    block::Int
end

"""
    CoincidenceWriter(path, attrs; block=2^20) -> CoincidenceWriter

Open `path` and create one extensible, chunked, shuffle+deflate dataset per LOR column, plus
the quantization scales and `attrs` (run metadata) as root attributes. `push_coincidence!`
buffers rows and flushes a block when `block` is reached; `close` flushes the remainder and
stamps `nrows`.
"""
function CoincidenceWriter(path::AbstractString, attrs::AbstractDict; block::Int = COINC_BLOCK)
    f = h5open(path, "w")
    proto = CoincidenceBuffer()
    dsets = Dict{String,HDF5.Dataset}()
    for (name, v) in coinc_columns(proto)
        dsets[name] = create_dataset(f, name, eltype(v), dataspace((0,); max_dims=(-1,));
                                     chunk=(block,), shuffle=true, deflate=4)
    end
    for (k, v) in attrs
        attributes(f)[String(k)] = v
    end
    attributes(f)["xyz_scale_mm"] = XYZ_SCALE_MM
    attributes(f)["e_scale_keV"]  = E_SCALE_KEV
    CoincidenceWriter(f, dsets, proto, 0, block)
end

function _flush!(w::CoincidenceWriter)
    m = length(w.buf)
    m == 0 && return
    newn = w.n + m
    for (name, v) in coinc_columns(w.buf)
        d = w.dsets[name]
        HDF5.set_extent_dims(d, (newn,))
        d[(w.n + 1):newn] = v
    end
    w.n = newn
    empty!(w.buf)
end

function push_coincidence!(w::CoincidenceWriter, args...)
    push_coincidence!(w.buf, args...)
    length(w.buf) >= w.block && _flush!(w)
    w
end

"Set a root attribute on the (open) LOR file — e.g. the final `nevents` once the stream ends."
set_lor_attr!(w::CoincidenceWriter, key::AbstractString, val) = (attributes(w.file)[key] = val; w)

function Base.close(w::CoincidenceWriter)::Int
    _flush!(w)
    attributes(w.file)["nrows"] = w.n
    close(w.file)
    w.n
end

"""
    foreach_coincidences_hdf5(f, path; batch=2^20) -> total

Stream a LOR HDF5 in row-slices, calling `f(buf::CoincidenceBuffer)` per slice (raw quantized
ints — decode with `decode_xyz`/`decode_e`). Returns the total row count.
"""
function foreach_coincidences_hdf5(f, path::AbstractString; batch::Int = 1 << 20)::Int
    h5open(path, "r") do h
        total = haskey(attributes(h), "nrows") ? Int(read(attributes(h)["nrows"])) : length(h["event"])
        lo = 1
        while lo <= total
            hi = min(lo + batch - 1, total); rng = lo:hi
            rd(name) = h[name][rng]
            b = CoincidenceBuffer(rd("event"), rd("truth"),
                rd("x1_mm"), rd("y1_mm"), rd("z1_mm"), rd("e1_keV"), rd("t1_ns"), rd("iz1"), rd("iphi1"),
                rd("x2_mm"), rd("y2_mm"), rd("z2_mm"), rd("e2_keV"), rd("t2_ns"), rd("iz2"), rd("iphi2"),
                rd("dt_ns"),
                rd("x0_mm"), rd("y0_mm"), rd("z0_mm"))
            f(b)
            lo = hi + 1
        end
        total
    end
end
