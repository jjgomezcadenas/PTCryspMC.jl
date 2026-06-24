# NIST XCOM photon cross-section loader + log-log interpolation, photon channels only.

"NIST XCOM photon cross sections; all columns mu/rho [cm^2/g]."
struct XCOMData
    E_MeV::Vector{Float64}
    coherent::Vector{Float64}
    incoherent::Vector{Float64}
    photoelectric::Vector{Float64}
    pair_nuclear::Vector{Float64}
    pair_electron::Vector{Float64}
    total_w_coh::Vector{Float64}
    total_no_coh::Vector{Float64}
end

"""
    load_xcom(path) -> XCOMData

Read a NIST XCOM table in the standard space-delimited format (2 header lines +
blank line + data rows of 8 numeric columns). Non-numeric lines are skipped.

Deferred (dev/status.md #2): harden before adding the crystal (CsI/BGO/LYSO) XCOM tables —
their in-range K-edge rows can carry a label that shifts columns, and empty input throws
obscurely. Validate the column count per row and error clearly on no data.
"""
function load_xcom(path::AbstractString)::XCOMData
    rows = Vector{Vector{Float64}}()
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            (s[1] >= '0' && s[1] <= '9') || continue
            fields = split(s)
            length(fields) >= 8 || continue
            push!(rows, parse.(Float64, fields[1:8]))
        end
    end
    m = reduce(hcat, rows)'
    XCOMData(m[:,1], m[:,2], m[:,3], m[:,4], m[:,5], m[:,6], m[:,7], m[:,8])
end

"Copy of the XCOM energy grid with duplicate K-edge energies nudged apart."
function _prepare_xcom_energy(xc::XCOMData)::Vector{Float64}
    E = copy(xc.E_MeV)
    @inbounds for i in 2:length(E)
        E[i] <= E[i-1] && (E[i] = E[i-1] * (1.0 + 1e-9))
    end
    E
end

"log.(fp), mapping non-positive values to -Inf (for interp_loglog_prelogged)."
prelog_data(fp::Vector{Float64})::Vector{Float64} = [f > 0.0 ? log(f) : -Inf for f in fp]

"""
    interp_loglog_prelogged(lx, log_xp, log_fp, fp, lo) -> Float64

Fast log-log interpolation using pre-computed log arrays and a known bracket
index `lo` (xp[lo] <= x < xp[lo+1]). One `exp` call; handles zero regions.
"""
function interp_loglog_prelogged(lx::Float64, log_xp::Vector{Float64},
                                 log_fp::Vector{Float64}, fp::Vector{Float64},
                                 lo::Int)::Float64
    hi = lo + 1
    @inbounds fp_lo = fp[lo]
    @inbounds fp_hi = fp[hi]
    fp_lo <= 0.0 && return 0.0                 # below-edge zero (also covers the both-zero case)
    fp_hi <= 0.0 && return fp_lo               # above-edge zero (e.g. pair below threshold) → flat
    @inbounds begin
        t = (lx - log_xp[lo]) / (log_xp[hi] - log_xp[lo])
        exp(log_fp[lo] + t * (log_fp[hi] - log_fp[lo]))
    end
end
