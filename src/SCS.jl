# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module SCS

import MathOptInterface
import Requires
import SCS_jll
import SparseArrays

function __init__()
    global indirect = SCS_jll.libscsindir
    global direct = SCS_jll.libscsdir
    Requires.@require(
        CUDA_jll = "e9e359dc-d701-5aa8-82ae-09bbf812ea83",
        begin
            import SCS_GPU_jll
            global gpuindirect = SCS_GPU_jll.libscsgpuindir
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
