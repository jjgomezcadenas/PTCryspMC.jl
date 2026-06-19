# HDF5 container for the singles stack: pack the binary part-files into one compressed HDF5,
# and stream it back for reading. Columnar typed datasets (chunked, shuffle + deflate-4 — the
# measured sweet spot: ~6× smaller than float CSV, fast typed reads), with the quantization
# scales and run metadata as root attributes (the self-describing `_meta`). All HDF5 calls are
# single-threaded here (the parallel chunks write plain binary parts) since the HDF5 C library
# is not thread-safe.

using HDF5

const H5_CHUNK   = 1 << 16     # dataset chunk length [rows]
const H5_DEFLATE = 4           # gzip level (deflate-9 buys ~1% for ~20× the write time)

"""
    pack_singles_hdf5(out, partpaths, rows, attrs) -> total

Read the binary part-files `partpaths` (chunk order, `rows[c]` rows each) and write `out`:
each singles column a chunked, shuffle+deflate dataset sized to the total, plus the
quantization scales and `attrs` (run metadata) as root attributes. Parts are streamed one at
a time (bounded memory) and deleted as consumed.
"""
function pack_singles_hdf5(out::AbstractString, partpaths::Vector{String}, rows::Vector{Int},
                           attrs::AbstractDict)::Int
    total = sum(rows)
    chunk = (max(1, min(total, H5_CHUNK)),)
    proto = SinglesBuffer()                                  # for the canonical (name, eltype) list
    h5open(out, "w") do f
        ds = Dict{String,Any}()
        for (name, v) in singles_columns(proto)
            ds[name] = create_dataset(f, name, eltype(v), (total,);
                                      chunk=chunk, shuffle=true, deflate=H5_DEFLATE)
        end
        attributes(f)["xyz_scale_mm"] = XYZ_SCALE_MM
        attributes(f)["e_scale_keV"]  = E_SCALE_KEV
        attributes(f)["nrows"]        = total
        for (k, v) in attrs
            attributes(f)[String(k)] = v
        end
        off = 0
        for c in eachindex(partpaths)
            n = rows[c]
            if n > 0
                buf = open(io -> read_part(io, n), partpaths[c], "r")
                rng = (off + 1):(off + n)
                for (name, v) in singles_columns(buf)
                    ds[name][rng] = v
                end
                off += n
            end
            rm(partpaths[c]; force=true)
        end
    end
    total
end

"Read a root attribute from a singles HDF5 (`default` if absent)."
singles_hdf5_attr(path::AbstractString, key::AbstractString, default=nothing) =
    h5open(h -> haskey(attributes(h), key) ? read(attributes(h)[key]) : default, path, "r")

"""
    foreach_singles_hdf5(f, path; batch=2^20) -> total

Stream a singles HDF5 in row-slices of at most `batch`, calling `f(buf::SinglesBuffer)` per
slice — bounded memory, for validation / downstream reads. Returns the total row count.
"""
function foreach_singles_hdf5(f, path::AbstractString; batch::Int = 1 << 20)::Int
    h5open(path, "r") do h
        total = haskey(attributes(h), "nrows") ? Int(read(attributes(h)["nrows"])) :
                                                 length(h["event_number"])
        lo = 1
        while lo <= total
            hi = min(lo + batch - 1, total)
            rng = lo:hi
            rd(name) = h[name][rng]
            buf = SinglesBuffer(rd("event_number"), rd("gamma"), rd("x_mm"), rd("y_mm"),
                                rd("z_mm"), rd("e_keV"), rd("iz"), rd("iphi"), rd("nblocks"),
                                rd("phantom_scatter"), rd("x0_mm"), rd("y0_mm"), rd("z0_mm"))
            f(buf)
            lo = hi + 1
        end
        total
    end
end
