struct DirectSolver <: LinearSolver end

scsint_t(::Type{DirectSolver}) = Int

function scs_set_default_settings(::Type{DirectSolver}, data::ScsData{Int})
    return ccall(
        (:scs_set_default_settings, direct),
        Cvoid,
        (Ptr{Cvoid},),
        data,
    )
end

function scs_init(
    ::Type{DirectSolver},
    data::ScsData{Int},
    cone::ScsCone{Int},
    info::ScsInfo{Int},
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

function scs_solve(
    ::Type{DirectSolver},
    p_work::Ptr{Cvoid},
    data::ScsData{Int},
    cone::ScsCone{Int},
    solution::ScsSolution,
    info::ScsInfo{Int},
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

function scs_finish(::Type{DirectSolver}, p_work::Ptr{Cvoid})
    return ccall((:scs_finish, direct), Cvoid, (Ptr{Cvoid},), p_work)
end

function scs_version()
    return unsafe_string(ccall((:scs_version, direct), Cstring, ()))
end
