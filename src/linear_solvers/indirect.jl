struct IndirectSolver <: LinearSolver end

scsint_t(::Type{IndirectSolver}) = Clonglong

function scs_set_default_settings(
    solver_t::Type{IndirectSolver},
    stgs::ScsSettings{I},
) where {I}
    @assert I == scsint_t(solver_t) == Clonglong
    return @ccall indirect.scs_set_default_settings(
        stgs::Ref{ScsSettings{I}},
    )::Cvoid
end

function scs_init(
    solver_t::Type{IndirectSolver},
    data::ScsData{I},
    cone::ScsCone{I},
    stgs::ScsSettings{I},
) where {I}
    @assert I == scsint_t(solver_t) == Clonglong
    return @ccall indirect.scs_init(
        data::Ref{ScsData{I}},
        cone::Ref{ScsCone{I}},
        stgs::Ref{ScsSettings{I}},
    )::Ptr{Cvoid}
end

function scs_solve(
    solver_t::Type{IndirectSolver},
    work::Ptr{Cvoid}, # ScsWork, unwrapped
    solution::ScsSolution,
    info::ScsInfo{I},
) where {I}
    @assert I == scsint_t(solver_t) == Clonglong
    return @ccall indirect.scs_solve(
        work::Ptr{Cvoid},
        solution::Ref{ScsSolution},
        info::Ref{ScsInfo},
    )::Clonglong
end

function scs_finish(::Type{IndirectSolver}, work::Ptr{Cvoid})
    return @ccall indirect.scs_finish(work::Ptr{Cvoid})::Cvoid
end
