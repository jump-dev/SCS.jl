# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

struct IndirectSolver <: LinearSolver end

is_available(::Type{IndirectSolver}) = true

scsint_t(::Type{IndirectSolver}) = Clonglong

function scs_set_default_settings(
    ::Type{IndirectSolver},
    stgs::ScsSettings{I},
) where {I<:Clonglong}
    return @ccall(
        libscsindir.scs_set_default_settings(stgs::Ref{ScsSettings{I}})::Cvoid,
    )
end

function scs_init(
    ::Type{IndirectSolver},
    data::ScsData{I},
    cone::ScsCone{I},
    stgs::ScsSettings{I},
) where {I<:Clonglong}
    return @ccall libscsindir.scs_init(
        data::Ref{ScsData{I}},
        cone::Ref{ScsCone{I}},
        stgs::Ref{ScsSettings{I}},
    )::Ptr{Cvoid}
end

function scs_update(
    ::Type{IndirectSolver},
    work::Ptr{Cvoid},
    b::Vector{Float64},
    c::Vector{Float64},
)
    return @ccall direct.scs_update(
        work::Ptr{Cvoid},
        b::Ref{Float64},
        c::Ref{Float64},
    )::Clonglong
end

function scs_solve(
    ::Type{IndirectSolver},
    work::Ptr{Cvoid},
    solution::ScsSolution,
    info::ScsInfo{I},
    warm_start::Integer,
) where {I<:Clonglong}
    return @ccall libscsindir.scs_solve(
        work::Ptr{Cvoid},
        solution::Ref{ScsSolution},
        info::Ref{ScsInfo{I}},
        warm_start::Clonglong,
    )::Clonglong
end

function scs_finish(::Type{IndirectSolver}, work::Ptr{Cvoid})
    return @ccall libscsindir.scs_finish(work::Ptr{Cvoid})::Cvoid
end

function scs_version(::Type{IndirectSolver})
    return unsafe_string(@ccall libscsindir.scs_version()::Cstring)
end
