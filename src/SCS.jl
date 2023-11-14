# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module SCS

import MathOptInterface as MOI
import SCS_jll
import SparseArrays

abstract type LinearSolver end

SCS.is_available(::Type{<:LinearSolver}) = false

include("c_wrapper.jl")
include("linear_solvers/direct.jl")
include("linear_solvers/indirect.jl")
include("MOI_wrapper/MOI_wrapper.jl")

# Code is contained in /ext/SCSSCS_GPU_jllExt
struct GpuIndirectSolver <: LinearSolver end

# Code is contained in /ext/SCSSCS_MKL_jllExt
struct MKLDirectSolver <: LinearSolver end

export scs_solve

end
