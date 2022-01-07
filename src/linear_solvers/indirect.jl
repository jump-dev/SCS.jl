struct IndirectSolver <: LinearSolver end

scsint_t(::Type{IndirectSolver}) = Clonglong

function scs_set_default_settings(
    ::Type{IndirectSolver},
    stgs::ScsSettings{I},
) where {I<:Clonglong}
    return @ccall(
        indirect.scs_set_default_settings(stgs::Ref{ScsSettings{I}})::Cvoid,
    )
end

function scs_init(
    ::Type{IndirectSolver},
    data::ScsData{I},
    cone::ScsCone{I},
    stgs::ScsSettings{I},
) where {I<:Clonglong}
    return @ccall indirect.scs_init(
        data::Ref{ScsData{I}},
        cone::Ref{ScsCone{I}},
        stgs::Ref{ScsSettings{I}},
    )::Ptr{Cvoid}
end

function scs_solve(
    ::Type{IndirectSolver},
    work::Ptr{Cvoid},
    solution::ScsSolution,
    info::ScsInfo{I},
) where {I<:Clonglong}
    return @ccall indirect.scs_solve(
        work::Ptr{Cvoid},
        solution::Ref{ScsSolution},
        info::Ref{ScsInfo{I}},
    )::Clonglong
end

function scs_finish(::Type{IndirectSolver}, work::Ptr{Cvoid})
    return @ccall indirect.scs_finish(work::Ptr{Cvoid})::Cvoid
end
