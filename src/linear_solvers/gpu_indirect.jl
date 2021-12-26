struct GpuIndirectSolver <: LinearSolver end

if haskey(ENV, "JULIA_SCS_LIBRARY_PATH")
    @isdefined(libscsgpuindir) && push!(available_solvers, GpuIndirectSolver)
else
    import SCS_GPU_jll
    const gpuindirect = SCS_GPU_jll.libscsgpuindir
    push!(available_solvers, GpuIndirectSolver)
end

scsint_t(::Type{GpuIndirectSolver}) = Cint

function scs_set_default_settings(
    solver_t::Type{GpuIndirectSolver},
    stgs::ScsSettings{I},
) where {I}
    @assert I == scsint_t(solver_t)
    return @ccall gpuindirect.scs_set_default_settings(
        stgs::Ref{ScsSettings{I}},
    )::Cvoid
end

function scs_init(
    solver_t::Type{GpuIndirectSolver},
    data::ScsData{I},
    cone::ScsCone{I},
    stgs::ScsSettings{I},
) where {I}
    @assert I == scsint_t(solver_t)
    return @ccall gpuindirect.scs_init(
        data::Ref{ScsData{I}},
        cone::Ref{ScsCone{I}},
        stgs::Ref{ScsSettings{I}},
    )::Ptr{Cvoid}
end

function scs_solve(
    solver_t::Type{GpuIndirectSolver},
    work::Ptr{Cvoid},  # ScsWork, unwrapped
    solution::ScsSolution,
    info::ScsInfo{I},
) where {I}
    @assert I == scsint_t(solver_t)
    return @ccall gpuindirect.scs_solve(
        work::Ptr{Cvoid},
        solution::Ref{ScsSolution},
        info::Ref{ScsInfo{I}},
    )::Cint
end

function scs_finish(::Type{GpuIndirectSolver}, work::Ptr{Cvoid})
    return @ccall gpuindirect.scs_finish(work::Ptr{Cvoid})::Cvoid
end
