struct GpuIndirectSolver <: LinearSolver end

scsint_t(::Type{GpuIndirectSolver}) = Cint

function scs_set_default_settings(
    ::Type{GpuIndirectSolver},
    stgs::ScsSettings{I},
) where {I<:Cint}
    return @ccall(
        gpuindirect.scs_set_default_settings(stgs::Ref{ScsSettings{I}})::Cvoid,
    )
end

function scs_init(
    ::Type{GpuIndirectSolver},
    data::ScsData{I},
    cone::ScsCone{I},
    stgs::ScsSettings{I},
) where {I<:Cint}
    return @ccall gpuindirect.scs_init(
        data::Ref{ScsData{I}},
        cone::Ref{ScsCone{I}},
        stgs::Ref{ScsSettings{I}},
    )::Ptr{Cvoid}
end

function scs_solve(
    ::Type{GpuIndirectSolver},
    work::Ptr{Cvoid},
    solution::ScsSolution,
    info::ScsInfo{I},
) where {I<:Cint}
    return @ccall gpuindirect.scs_solve(
        work::Ptr{Cvoid},
        solution::Ref{ScsSolution},
        info::Ref{ScsInfo{I}},
    )::Cint
end

function scs_finish(::Type{GpuIndirectSolver}, work::Ptr{Cvoid})
    return @ccall gpuindirect.scs_finish(work::Ptr{Cvoid})::Cvoid
end
