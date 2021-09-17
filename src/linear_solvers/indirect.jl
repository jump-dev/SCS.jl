struct IndirectSolver <: LinearSolver end

scsint_t(::Type{IndirectSolver}) = Int

function scs_set_default_settings(
    ::Type{IndirectSolver},
    stgs::ScsSettings{Int},
)
    return ccall(
        (:scs_set_default_settings, indirect),
        Cvoid,
        (Ref{ScsSettings{Int}},),
        stgs,
    )
end

function scs_init(
    ::Type{IndirectSolver},
    data::ScsData{Int},
    cone::ScsCone{Int},
    stgs::ScsSettings{Int},
)
    return ccall(
        (:scs_init, indirect),
        Ptr{Cvoid}, # returns ptr to unwrapped ScsWork
        (Ref{ScsData{Int}}, Ref{ScsCone{Int}}, Ref{ScsSettings{Int}}),
        data,
        cone,
        stgs,
    )
end

function scs_solve(
    ::Type{IndirectSolver},
    work::Ptr{Cvoid}, # ScsWork, unwrapped
    solution::ScsSolution,
    info::ScsInfo{Int},
)
    return ccall(
        (:scs_solve, indirect),
        Int,
        (Ptr{Cvoid}, Ref{ScsSolution}, Ref{ScsInfo{Int}}),
        work,
        solution,
        info,
    )
end

function scs_finish(::Type{IndirectSolver}, work::Ptr{Cvoid})
    return ccall((:scs_finish, indirect), Cvoid, (Ptr{Cvoid},), work)
end
