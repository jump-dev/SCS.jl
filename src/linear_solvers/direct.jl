struct DirectSolver <: LinearSolver end

scsint_t(::Type{DirectSolver}) = Clonglong

function scs_set_default_settings(
    solver_t::Type{DirectSolver},
    stgs::ScsSettings{I},
) where {I}
    @assert I == scsint_t(solver_t)
    return @ccall direct.scs_set_default_settings(
        stgs::Ref{ScsSettings{I}},
    )::Cvoid
end

function scs_init(
    solver_t::Type{DirectSolver},
    data::ScsData{I},
    cone::ScsCone{I},
    stgs::ScsSettings{I},
) where {I}
    @assert I == scsint_t(solver_t)
    return @ccall direct.scs_init(
        data::Ref{ScsData{I}},
        cone::Ref{ScsCone{I}},
        stgs::Ref{ScsSettings{I}},
    )::Ptr{Cvoid}
end

function scs_solve(
    solver_t::Type{DirectSolver},
    work::Ptr{Cvoid}, # ScsWork, unwrapped
    solution::ScsSolution,
    info::ScsInfo{I},
) where {I}
    @assert I == scsint_t(solver_t)
    return @ccall direct.scs_solve(
        work::Ptr{Cvoid},
        solution::Ref{ScsSolution},
        info::Ref{ScsInfo{I}},
    )::Clonglong
end

function scs_finish(::Type{DirectSolver}, work::Ptr{Cvoid})
    return @ccall direct.scs_finish(work::Ptr{Cvoid})::Cvoid
end

function scs_version()
    return unsafe_string(@ccall direct.scs_version()::Cstring)
end
