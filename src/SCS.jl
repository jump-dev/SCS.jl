# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module SCS

import MathOptInterface
import Requires
import SCS_jll
import SparseArrays

global indirect = SCS_jll.libscsindir
global direct = SCS_jll.libscsdir

include("c_wrapper.jl")
include("linear_solvers/direct.jl")
include("linear_solvers/indirect.jl")
include("MOI_wrapper/MOI_wrapper.jl")

const available_solvers = [DirectSolver, IndirectSolver]

struct GpuIndirectSolver <: LinearSolver end
struct MKLDirectSolver <: LinearSolver end

export scs_solve

end
