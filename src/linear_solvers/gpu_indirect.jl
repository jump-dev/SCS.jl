struct GpuIndirectSolver <: LinearSolver end

if haskey(ENV, "JULIA_SCS_LIBRARY_PATH")
    @isdefined(libscsgpuindir) && push!(available_solvers, GpuIndirectSolver)
else
    import SCS_GPU_jll
    const gpuindirect = SCS_GPU_jll.libscsgpuindir
    push!(available_solvers, GpuIndirectSolver)
end

scsint_t(::Type{GpuIndirectSolver}) = Cint

function SCS_set_default_settings(
    ::Type{GpuIndirectSolver},
    data::SCSData{Cint},
)
    return ccall(
        (:scs_set_default_settings, gpuindirect),
        Cvoid,
        (Ptr{Cvoid},),
        data,
    )
end

function SCS_init(
    ::Type{GpuIndirectSolver},
    data::SCSData{Cint},
    cone::SCSCone{Cint},
    info::SCSInfo{Cint},
)
    return ccall(
        (:scs_init, gpuindirect),
        Ptr{Cvoid},
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        data,
        cone,
        info,
    )
end

function SCS_solve(
    ::Type{GpuIndirectSolver},
    p_work::Ptr{Cvoid},
    data::SCSData{Cint},
    cone::SCSCone{Cint},
    solution::SCSSolution,
    info::SCSInfo{Cint},
)
    return ccall(
        (:scs_solve, gpuindirect),
        Cint,
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        p_work,
        data,
        cone,
        solution,
        info,
    )
end

function SCS_finish(::Type{GpuIndirectSolver}, p_work::Ptr{Cvoid})
    return ccall((:scs_finish, gpuindirect), Cvoid, (Ptr{Cvoid},), p_work)
end
