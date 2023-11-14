# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module SCS

import MathOptInterface as MOI
import SCS_jll
import SparseArrays

abstract type LinearSolver end

is_available(::Type{<:LinearSolver}) = false

include("c_wrapper.jl")
include("linear_solvers/direct.jl")
include("linear_solvers/indirect.jl")
include("MOI_wrapper/MOI_wrapper.jl")

# Code is contained in /ext/SCSSCS_GPU_jllExt
struct GpuIndirectSolver <: LinearSolver end

# Code is contained in /ext/SCSSCS_MKL_jllExt
struct MKLDirectSolver <: LinearSolver end

@static if !isdefined(Base, :get_extension)
    import Requires
    function __init__()
        Requires.@require(
            SCS_GPU_jll = "af6e375f-46ec-5fa0-b791-491b0dfa44a4",
            include("../ext/SCSSCS_GPU_jllExt.jl"),
        )
        Requires.@require(
            SCS_MKL_jll = "3f2553a9-4106-52be-b7dd-865123654657",
            include("../ext/SCSSCS_MKL_jllExt.jl"),
        )
        return
    end
end

export scs_solve

end
