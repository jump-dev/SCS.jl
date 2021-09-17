struct DirectSolver <: LinearSolver end

scsint_t(::Type{DirectSolver}) = Int

function scs_set_default_settings(::Type{DirectSolver}, stgs::ScsSettings{Int})
    return ccall(
        (:scs_set_default_settings, direct),
        Cvoid,
        (Ref{ScsSettings{Int}},),
        stgs,
    )
end

function scs_init(
    ::Type{DirectSolver},
    data::ScsData{Int},
    cone::ScsCone{Int},
    stgs::ScsSettings{Int},
)
    return ccall(
        (:scs_init, direct),
        Ptr{Cvoid}, # returns ptr to unwrapped ScsWork
        (Ref{ScsData{Int}}, Ref{ScsCone{Int}}, Ref{ScsSettings{Int}}),
        data,
        cone,
        stgs,
    )
end

function scs_solve(
    ::Type{DirectSolver},
    work::Ptr{Cvoid}, # ScsWork, unwrapped
    solution::ScsSolution,
    info::ScsInfo{Int},
)
    return ccall(
        (:scs_solve, direct),
        Int,
        (Ptr{Cvoid}, Ref{ScsSolution}, Ref{ScsInfo{Int}}),
        work,
        solution,
        info,
    )
end

function scs_finish(::Type{DirectSolver}, work::Ptr{Cvoid})
    return ccall((:scs_finish, direct), Cvoid, (Ptr{Cvoid},), work)
end

function scs_version()
    return unsafe_string(ccall((:scs_version, direct), Cstring, ()))
end
