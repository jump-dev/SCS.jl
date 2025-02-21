# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

struct DirectSolver <: LinearSolver end

is_available(::Type{DirectSolver}) = true

scsint_t(::Type{DirectSolver}) = Clonglong

function scs_set_default_settings(
    ::Type{DirectSolver},
    stgs::ScsSettings{I},
) where {I<:Clonglong}
    return @ccall(
        libscsdir.scs_set_default_settings(stgs::Ref{ScsSettings{I}})::Cvoid,
    )
end

function scs_init(
    ::Type{DirectSolver},
    data::ScsData{I},
    cone::ScsCone{I},
    stgs::ScsSettings{I},
) where {I<:Clonglong}
    return @ccall libscsdir.scs_init(
        data::Ref{ScsData{I}},
        cone::Ref{ScsCone{I}},
        stgs::Ref{ScsSettings{I}},
    )::Ptr{Cvoid}
end

function scs_update(
    ::Type{DirectSolver},
    work::Ptr{Cvoid},
    b::Vector{Float64},
    c::Vector{Float64},
)
    return @ccall libscsdir.scs_update(
        work::Ptr{Cvoid},
        b::Ref{Float64},
        c::Ref{Float64},
    )::Clonglong
end

function scs_solve(
    ::Type{DirectSolver},
    work::Ptr{Cvoid},
    solution::ScsSolution,
    info::ScsInfo{I},
    warm_start::Integer,
) where {I<:Clonglong}
    return @ccall libscsdir.scs_solve(
        work::Ptr{Cvoid},
        solution::Ref{ScsSolution},
        info::Ref{ScsInfo{I}},
        warm_start::Clonglong,
    )::Clonglong
end

function scs_finish(::Type{DirectSolver}, work::Ptr{Cvoid})
    return @ccall libscsdir.scs_finish(work::Ptr{Cvoid})::Cvoid
end

function scs_version(::Type{DirectSolver})
    return unsafe_string(@ccall libscsdir.scs_version()::Cstring)
end

scs_version() = scs_version(DirectSolver)
