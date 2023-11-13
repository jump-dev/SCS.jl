# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module SCSSCS_GPU_jllExt

import SCS
import SCS_GPU_jll

global gpuindirect = SCS_GPU_jll.libscsgpuindir

SCS.is_available(::Type{SCS.GpuIndirectSolver}) = true

SCS.scsint_t(::Type{SCS.GpuIndirectSolver}) = Cint

function SCS.scs_set_default_settings(
    ::Type{SCS.GpuIndirectSolver},
    stgs::SCS.ScsSettings{I},
) where {I<:Cint}
    return @ccall gpuindirect.scs_set_default_settings(
        stgs::Ref{SCS.ScsSettings{I}},
    )::Cvoid,
end

function SCS.scs_init(
    ::Type{SCS.GpuIndirectSolver},
    data::SCS.ScsData{I},
    cone::SCS.ScsCone{I},
    stgs::SCS.ScsSettings{I},
) where {I<:Cint}
    return @ccall gpuindirect.scs_init(
        data::Ref{SCS.ScsData{I}},
        cone::Ref{SCS.ScsCone{I}},
        stgs::Ref{SCS.ScsSettings{I}},
    )::Ptr{Cvoid}
end

function SCS.scs_update(
    ::Type{SCS.GpuIndirectSolver},
    work::Ptr{Cvoid},
    b::Vector{Float64},
    c::Vector{Float64},
)
    return @ccall direct.scs_update(
        work::Ptr{Cvoid},
        b::Ref{Float64},
        c::Ref{Float64},
    )::Cint
end

function SCS.scs_solve(
    ::Type{SCS.GpuIndirectSolver},
    work::Ptr{Cvoid},
    solution::SCS.ScsSolution,
    info::SCS.ScsInfo{I},
    warm_start::Integer,
) where {I<:Cint}
    return @ccall gpuindirect.scs_solve(
        work::Ptr{Cvoid},
        solution::Ref{SCS.ScsSolution},
        info::Ref{SCS.ScsInfo{I}},
        warm_start::Cint,
    )::Cint
end

function SCS.scs_finish(::Type{SCS.GpuIndirectSolver}, work::Ptr{Cvoid})
    return @ccall gpuindirect.scs_finish(work::Ptr{Cvoid})::Cvoid
end

function SCS.scs_version(::Type{SCS.GpuIndirectSolver})
    return unsafe_string(@ccall gpuindirect.scs_version()::Cstring)
end

end  # module
