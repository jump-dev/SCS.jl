module SCS

import Requires
import SCS_jll
import SparseArrays

const indirect = SCS_jll.libscsindir
const direct = SCS_jll.libscsdir

function __init__()
    Requires.@require(
        CUDA_jll = "e9e359dc-d701-5aa8-82ae-09bbf812ea83",
        begin
            import SCS_GPU_jll
            const gpuindirect = SCS_GPU_jll.libscsgpuindir
            push!(available_solvers, GpuIndirectSolver)
        end
    )
    return
end

include("c_wrapper.jl")
include("linear_solvers/direct.jl")
include("linear_solvers/indirect.jl")
include("linear_solvers/gpu_indirect.jl")
include("MOI_wrapper/MOI_wrapper.jl")

const available_solvers = [DirectSolver, IndirectSolver]

export scs_solve

end
