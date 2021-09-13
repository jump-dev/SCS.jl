struct IndirectSolver <: LinearSolver end

scsint_t(::Type{IndirectSolver}) = Int

function SCS_set_default_settings(::Type{IndirectSolver}, data::SCSData{Int})
    return ccall(
        (:scs_set_default_settings, indirect),
        Cvoid,
        (Ptr{Cvoid},),
        data,
    )
end

function SCS_init(
    ::Type{IndirectSolver},
    data::SCSData{Int},
    cone::SCSCone{Int},
    info::SCSInfo{Int},
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

function SCS_solve(
    ::Type{IndirectSolver},
    p_work::Ptr{Cvoid},
    data::SCSData{Int},
    cone::SCSCone{Int},
    solution::SCSSolution,
    info::SCSInfo{Int},
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

function SCS_finish(::Type{IndirectSolver}, p_work::Ptr{Cvoid})
    return ccall((:scs_finish, indirect), Cvoid, (Ptr{Cvoid},), p_work)
end
