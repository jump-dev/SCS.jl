struct IndirectSolver <: LinearSolver end

scsint_t(::Type{IndirectSolver}) = Int

function scs_set_default_settings(::Type{IndirectSolver}, data::ScsData{Int})
    return ccall(
        (:scs_set_default_settings, indirect),
        Cvoid,
        (Ptr{Cvoid},),
        data,
    )
end

function scs_init(
    ::Type{IndirectSolver},
    data::ScsData{Int},
    cone::ScsCone{Int},
    info::ScsInfo{Int},
)
    return ccall(
        (:scs_init, indirect),
        Ptr{Cvoid},
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        data,
        cone,
        info,
    )
end

function scs_solve(
    ::Type{IndirectSolver},
    p_work::Ptr{Cvoid},
    data::ScsData{Int},
    cone::ScsCone{Int},
    solution::ScsSolution,
    info::ScsInfo{Int},
)
    return ccall(
        (:scs_solve, indirect),
        Int,
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        p_work,
        data,
        cone,
        solution,
        info,
    )
end

function scs_finish(::Type{IndirectSolver}, p_work::Ptr{Cvoid})
    return ccall((:scs_finish, indirect), Cvoid, (Ptr{Cvoid},), p_work)
end
