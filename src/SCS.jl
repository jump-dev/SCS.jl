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

export scs_solve

import PrecompileTools

PrecompileTools.@setup_workload begin
    PrecompileTools.@compile_workload begin
        model = MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            MOI.instantiate(SCS.Optimizer; with_bridge_type = Float64),
        )
        MOI.set(model, MOI.Silent(), true)
        x = MOI.add_variables(model, 3)
        MOI.supports(model, MOI.VariableName(), typeof(x[1]))
        MOI.set(model, MOI.VariableName(), x[1], "x1")
        MOI.set(model, MOI.VariablePrimalStart(), x[1], 0.0)
        sets = (MOI.GreaterThan(0.0), MOI.LessThan(0.0), MOI.GreaterThan(0.0))
        for (i, set) in enumerate(sets)
            MOI.add_constrained_variable(model, set)
            MOI.supports_constraint(model, typeof(x[i]), typeof(set))
            MOI.add_constraint(model, x[i], set)
            MOI.supports_constraint(model, typeof(1.0 * x[i]), typeof(set))
            c = MOI.add_constraint(model, 1.0 * x[i] + 0.0, set)
            MOI.supports(model, MOI.ConstraintName(), typeof(c))
            MOI.set(model, MOI.ConstraintName(), c, "c1")
        end
        f = 1.0 * x[1] + x[2] + x[3]
        MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(model, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        g_vov = MOI.VectorOfVariables(x)
        g_vaf = MOI.Utilities.vectorize(1.0 .* x)
        for set in (MOI.ExponentialCone(), MOI.SecondOrderCone(3))
            MOI.add_constraint(model, g_vov, set)
            MOI.add_constraint(model, g_vaf, set)
        end
        MOI.supports(model, MOI.ObjectiveFunction{typeof(f)}())
        MOI.set(model, MOI.ObjectiveFunction{typeof(f)}(), f)
        MOI.optimize!(model)
        MOI.get(model, MOI.TerminationStatus())
        MOI.get(model, MOI.PrimalStatus())
        MOI.get(model, MOI.DualStatus())
        MOI.get(model, MOI.VariablePrimal(), x)
    end
end

end
