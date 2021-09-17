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
    stgs::ScsSettings{Int32},
)
    return ccall(
        (:scs_set_default_settings, gpuindirect),
        Cvoid,
        (Ref{ScsSettings{Int32}},),
        stgs,
    )
end

function scs_init(
    ::Type{GpuIndirectSolver},
    data::ScsData{Int32},
    cone::ScsCone{Int32},
    stgs::ScsSettings{Int32},
)
    return ccall(
        (:scs_init, gpuindirect),
        Ptr{Cvoid}, # returns ptr to unwrapped ScsWork
        (Ref{ScsData{Int32}}, Ref{ScsCone{Int32}}, Ref{ScsSettings{Int32}}),
        data,
        cone,
        stgs,
    )
end

function scs_solve(
    ::Type{GpuIndirectSolver},
    work::Ptr{Cvoid},  # ScsWork, unwrapped
    solution::ScsSolution,
    info::ScsInfo{Int32},
)
    return ccall(
        (:scs_solve, gpuindirect),
        Int32,
        (Ptr{Cvoid}, Ref{ScsSolution}, Ref{ScsInfo{Int32}}),
        work,
        solution,
        info,
    )
end

function scs_finish(::Type{GpuIndirectSolver}, work::Ptr{Cvoid})
    return ccall((:scs_finish, gpuindirect), Cvoid, (Ptr{Cvoid},), work)
end
