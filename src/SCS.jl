module SCS

using Libdl
import Requires
import SparseArrays

if haskey(ENV, "JULIA_SCS_LIBRARY_PATH")
    if !isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
        error(
            "SCS not properly installed. Please run `Pkg.build(\"SCS\")` and " *
            "restart julia",
        )
    end
    include(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    function __init__()
        version = VersionNumber(scs_version())
        if version < v"3.0.0"
            error(
                "Current SCS version installed is $version, but we require " *
                "version 3.0.*",
            )
        end
        return
    end
else
    import SCS_jll
    const indirect = SCS_jll.libscsindir
    const direct = SCS_jll.libscsdir
    function __init__()
        Requires.@require(
            CUDA_jll = "e9e359dc-d701-5aa8-82ae-09bbf812ea83",
            if haskey(ENV, "JULIA_SCS_LIBRARY_PATH")
                @isdefined(libscsgpuindir) &&
                    push!(available_solvers, GpuIndirectSolver)
            else
                import SCS_GPU_jll
                const gpuindirect = SCS_GPU_jll.libscsgpuindir
                push!(available_solvers, GpuIndirectSolver)
            end
        )
        return
    end
end

include("c_wrapper.jl")
include("linear_solvers/direct.jl")
include("linear_solvers/indirect.jl")
include("linear_solvers/gpu_indirect.jl")
include("MOI_wrapper/MOI_wrapper.jl")

const available_solvers = [DirectSolver, IndirectSolver]

export scs_solve

end
