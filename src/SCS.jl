module SCS

using Libdl
import Requires
import SparseArrays

if haskey(ENV, "JULIA_SCS_LIBRARY_PATH") || VERSION < v"1.3"
    if !isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
        error(
            "SCS not properly installed. Please run `Pkg.build(\"SCS\")` and " *
            "restart julia",
        )
    end
    include(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    function __init__()
        version = VersionNumber(SCS_version())
        if version < v"2.1.0"
            error(
                "Current SCS version installed is $version, but we require " *
                "version 2.1.*",
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
            include("linear_solvers/gpu_indirect.jl"),
        )
        return
    end
end

include("c_wrapper.jl")
include("linear_solvers/direct.jl")
include("linear_solvers/indirect.jl")
include("MOI_wrapper/MOI_wrapper.jl")

const available_solvers = [DirectSolver, IndirectSolver]

export SCS_solve

end
