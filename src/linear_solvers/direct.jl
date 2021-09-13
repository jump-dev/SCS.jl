struct DirectSolver <: LinearSolver end

scsint_t(::Type{DirectSolver}) = Int

function SCS_set_default_settings(::Type{DirectSolver}, data::SCSData{Int})
    return ccall(
        (:scs_set_default_settings, direct),
        Cvoid,
        (Ptr{Cvoid},),
        data,
    )
end

function SCS_init(
    ::Type{DirectSolver},
    data::SCSData{Int},
    cone::SCSCone{Int},
    info::SCSInfo{Int},
)
    return ccall(
        (:scs_init, direct),
        Ptr{Cvoid},
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        data,
        cone,
        info,
    )
end

function SCS_solve(
    ::Type{DirectSolver},
    p_work::Ptr{Cvoid},
    data::SCSData{Int},
    cone::SCSCone{Int},
    solution::SCSSolution,
    info::SCSInfo{Int},
)
    return ccall(
        (:scs_solve, direct),
        Int,
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        p_work,
        data,
        cone,
        solution,
        info,
    )
end

function SCS_finish(::Type{DirectSolver}, p_work::Ptr{Cvoid})
    return ccall((:scs_finish, direct), Cvoid, (Ptr{Cvoid},), p_work)
end

function SCS_version()
    return unsafe_string(ccall((:scs_version, direct), Cstring, ()))
end
