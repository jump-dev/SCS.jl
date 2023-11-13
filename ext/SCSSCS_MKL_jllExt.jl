# Copyright (c) 2022: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module SCSSCS_MKL_jllExt

import SCS
import SCS_MKL_jll

function __init__()
    global mkldirect = SCS_MKL_jll.libscsmkl
    push!(SCS.available_solvers, SCS.MKLDirectSolver)
    return
end

SCS.scsint_t(::Type{SCS.MKLDirectSolver}) = Clonglong

function SCS.scs_set_default_settings(
    ::Type{SCS.MKLDirectSolver},
    stgs::SCS.ScsSettings{I},
) where {I<:Clonglong}
    return @ccall(
        mkldirect.scs_set_default_settings(stgs::Ref{SCS.ScsSettings{I}})::Cvoid,
    )
end

function SCS.scs_init(
    ::Type{SCS.MKLDirectSolver},
    data::SCS.ScsData{I},
    cone::SCS.ScsCone{I},
    stgs::SCS.ScsSettings{I},
) where {I<:Clonglong}
    return @ccall mkldirect.scs_init(
        data::Ref{SCS.ScsData{I}},
        cone::Ref{SCS.ScsCone{I}},
        stgs::Ref{SCS.ScsSettings{I}},
    )::Ptr{Cvoid}
end

function SCS.scs_update(
    ::Type{SCS.MKLDirectSolver},
    work::Ptr{Cvoid},
    b::Vector{Float64},
    c::Vector{Float64},
)
    return @ccall mkldirect.scs_update(
        work::Ptr{Cvoid},
        b::Ref{Float64},
        c::Ref{Float64},
    )::Clonglong
end

function SCS.scs_solve(
    ::Type{SCS.MKLDirectSolver},
    work::Ptr{Cvoid},
    solution::SCS.ScsSolution,
    info::SCS.ScsInfo{I},
    warm_start::Integer,
) where {I<:Clonglong}
    return @ccall mkldirect.scs_solve(
        work::Ptr{Cvoid},
        solution::Ref{SCS.ScsSolution},
        info::Ref{SCS.ScsInfo{I}},
        warm_start::Clonglong,
    )::Clonglong
end

function SCS.scs_finish(::Type{SCS.MKLDirectSolver}, work::Ptr{Cvoid})
    return @ccall mkldirect.scs_finish(work::Ptr{Cvoid})::Cvoid
end

function SCS.scs_version(::Type{SCS.MKLDirectSolver})
    return unsafe_string(@ccall mkldirect.scs_version()::Cstring)
end

end  # module
