struct GpuIndirectSolver <: LinearSolver end

if haskey(ENV, "JULIA_SCS_LIBRARY_PATH")
    @isdefined(libscsgpuindir) && push!(available_solvers, GpuIndirectSolver)
else
    import SCS_GPU_jll
    const gpuindirect = SCS_GPU_jll.libscsgpuindir
    push!(available_solvers, GpuIndirectSolver)
end

scsint_t(::Type{GpuIndirectSolver}) = Int32

function scs_set_default_settings(
    ::Type{GpuIndirectSolver},
    data::ScsData{Int32},
)
    return ccall(
        (:scs_set_default_settings, gpuindirect),
        Cvoid,
        (Ptr{Cvoid},),
        data,
    )
end

function scs_init(
    ::Type{GpuIndirectSolver},
    data::ScsData{Int32},
    cone::ScsCone{Int32},
    info::ScsInfo{Int32},
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

function scs_solve(
    ::Type{GpuIndirectSolver},
    p_work::Ptr{Cvoid},
    data::ScsData{Int32},
    cone::ScsCone{Int32},
    solution::ScsSolution,
    info::ScsInfo{Int32},
)
    return ccall(
        (:scs_solve, gpuindirect),
        Int32,
        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        p_work,
        data,
        cone,
        solution,
        info,
    )
end

function scs_finish(::Type{GpuIndirectSolver}, p_work::Ptr{Cvoid})
    return ccall((:scs_finish, gpuindirect), Cvoid, (Ptr{Cvoid},), p_work)
end
