# Copyright (c) 2022: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

struct MKLDirectSolver <: LinearSolver end

scsint_t(::Type{MKLDirectSolver}) = Clonglong

function scs_set_default_settings(
    ::Type{MKLDirectSolver},
    stgs::ScsSettings{I},
) where {I<:Clonglong}
    return @ccall(
        mkldirect.scs_set_default_settings(stgs::Ref{ScsSettings{I}})::Cvoid,
    )
end

function scs_init(
    ::Type{MKLDirectSolver},
    data::ScsData{I},
    cone::ScsCone{I},
    stgs::ScsSettings{I},
) where {I<:Clonglong}
    return @ccall mkldirect.scs_init(
        data::Ref{ScsData{I}},
        cone::Ref{ScsCone{I}},
        stgs::Ref{ScsSettings{I}},
    )::Ptr{Cvoid}
end

function scs_update(
    ::Type{MKLDirectSolver},
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

function scs_solve(
    ::Type{MKLDirectSolver},
    work::Ptr{Cvoid},
    solution::ScsSolution,
    info::ScsInfo{I},
    warm_start::Integer,
) where {I<:Clonglong}
    return @ccall mkldirect.scs_solve(
        work::Ptr{Cvoid},
        solution::Ref{ScsSolution},
        info::Ref{ScsInfo{I}},
        warm_start::Clonglong,
    )::Clonglong
end

function scs_finish(::Type{MKLDirectSolver}, work::Ptr{Cvoid})
    return @ccall mkldirect.scs_finish(work::Ptr{Cvoid})::Cvoid
end

function scs_version(::Type{MKLDirectSolver})
    return unsafe_string(@ccall mkldirect.scs_version()::Cstring)
end
