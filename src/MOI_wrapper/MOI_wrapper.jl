# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

include("scaled_psd_cone_bridge.jl")

MOI.Utilities.@product_of_sets(
    _Cones,
    MOI.Zeros,
    MOI.Nonnegatives,
    MOI.SecondOrderCone,
    ScaledPSDCone,
    MOI.ExponentialCone,
    MOI.DualExponentialCone,
    MOI.PowerCone{T},
    MOI.DualPowerCone{T},
)

struct _SetConstants{T}
    b::Vector{T}
    power_coefficients::Dict{Int,T}
    _SetConstants{T}() where {T} = new{T}(T[], Dict{Int,T}())
end

function Base.empty!(x::_SetConstants)
    empty!(x.b)
    empty!(x.power_coefficients)
    return x
end

Base.resize!(x::_SetConstants, n) = resize!(x.b, n)

function MOI.Utilities.load_constants(x::_SetConstants, offset, f)
    MOI.Utilities.load_constants(x.b, offset, f)
    return
end

function MOI.Utilities.load_constants(
    x::_SetConstants{T},
    offset,
    set::Union{MOI.PowerCone{T},MOI.DualPowerCone{T}},
) where {T}
    x.power_coefficients[offset+1] = set.exponent
    return
end

function MOI.Utilities.function_constants(x::_SetConstants, rows)
    return MOI.Utilities.function_constants(x.b, rows)
end

function MOI.Utilities.set_from_constants(x::_SetConstants, S, rows)
    return MOI.Utilities.set_from_constants(x.b, S, rows)
end

function MOI.Utilities.set_from_constants(
    x::_SetConstants{T},
    ::Type{S},
    rows,
) where {T,S<:Union{MOI.PowerCone{T},MOI.DualPowerCone{T}}}
    @assert length(rows) == 3
    return S(x.power_coefficients[first(rows)])
end

const OptimizerCache{T} = MOI.Utilities.GenericModel{
    Cdouble,
    MOI.Utilities.ObjectiveContainer{Cdouble},
    MOI.Utilities.VariablesContainer{Cdouble},
    MOI.Utilities.MatrixOfConstraints{
        Cdouble,
        MOI.Utilities.MutableSparseMatrixCSC{
            Cdouble,
            T,
            MOI.Utilities.ZeroBasedIndexing,
        },
        _SetConstants{Cdouble},
        _Cones{Cdouble},
    },
}

function _to_sparse(
    ::Type{T},
    A::MOI.Utilities.MutableSparseMatrixCSC{
        Cdouble,
        T,
        MOI.Utilities.ZeroBasedIndexing,
    },
) where {T}
    return -A.nzval, A.rowval, A.colptr
end

mutable struct MOISolution
    primal::Vector{Float64}
    dual::Vector{Float64}
    slack::Vector{Float64}
    ret_val::Int
    raw_status::String
    objective_value::Float64
    dual_objective_value::Float64
    objective_constant::Float64
    solve_time_sec::Float64
    iterations::Int
end

function MOISolution()
    return MOISolution(
        Float64[],
        Float64[],
        Float64[],
        0, # SCS_UNFINISHED
        "",
        NaN,
        NaN,
        NaN,
        0.0,
        0,
    )
end

"""
    Optimizer()

Create a new SCS optimizer.
"""
mutable struct Optimizer <: MOI.AbstractOptimizer
    cones::Union{Nothing,_Cones{Cdouble}}
    sol::MOISolution
    silent::Bool
    options::Dict{Symbol,Any}

    Optimizer() = new(nothing, MOISolution(), false, Dict{Symbol,Any}())
end

function MOI.get(::Optimizer, ::MOI.Bridges.ListOfNonstandardBridges)
    return [ScaledPSDConeBridge{Cdouble}]
end

MOI.get(::Optimizer, ::MOI.SolverName) = "SCS"

MOI.is_empty(optimizer::Optimizer) = optimizer.cones === nothing

function MOI.empty!(optimizer::Optimizer)
    optimizer.cones = nothing
    optimizer.sol = MOISolution()
    return
end

###
### MOI.RawOptimizerAttribute
###

function MOI.set(optimizer::Optimizer, param::MOI.RawOptimizerAttribute, value)
    return optimizer.options[Symbol(param.name)] = value
end

function MOI.get(optimizer::Optimizer, param::MOI.RawOptimizerAttribute)
    return optimizer.options[Symbol(param.name)]
end

###
### MOI.Silent
###

MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    optimizer.silent = value
    return
end

MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent

###
### MOI.TimeLimitSec
###

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(optimizer::Optimizer, ::MOI.TimeLimitSec, time_limit::Real)
    optimizer.options[:time_limit_secs] = convert(Float64, time_limit)
    return
end

function MOI.set(optimizer::Optimizer, ::MOI.TimeLimitSec, ::Nothing)
    delete!(optimizer.options, :time_limit_secs)
    return
end

function MOI.get(optimizer::Optimizer, ::MOI.TimeLimitSec)
    value = get(optimizer.options, :time_limit_secs, nothing)
    return value::Union{Float64,Nothing}
end

###
### MOI.AbstractModelAttribute
###

function MOI.supports(
    ::Optimizer,
    ::Union{
        MOI.ObjectiveSense,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
        MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Float64}},
    },
)
    return true
end

###
### MOI.AbstractVariableAttribute
###

function MOI.supports(
    ::Optimizer,
    ::MOI.VariablePrimalStart,
    ::Type{MOI.VariableIndex},
)
    return true
end

###
### MOI.AbstractConstraintAttribute
###

function MOI.supports(
    ::Optimizer,
    ::Union{MOI.ConstraintPrimalStart,MOI.ConstraintDualStart},
    ::Type{<:MOI.ConstraintIndex},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VectorAffineFunction{Cdouble}},
    ::Type{
        <:Union{
            MOI.Zeros,
            MOI.Nonnegatives,
            MOI.SecondOrderCone,
            ScaledPSDCone,
            MOI.ExponentialCone,
            MOI.DualExponentialCone,
            MOI.PowerCone{Cdouble},
            MOI.DualPowerCone{Cdouble},
        },
    },
)
    return true
end

function _map_sets(f, ::Type{T}, sets, ::Type{S}) where {T,S}
    F = MOI.VectorAffineFunction{Cdouble}
    cis = MOI.get(sets, MOI.ListOfConstraintIndices{F,S}())
    return T[f(MOI.get(sets, MOI.ConstraintSet(), ci)) for ci in cis]
end

function MOI.optimize!(
    dest::Optimizer,
    src::MOI.Utilities.UniversalFallback{OptimizerCache{T}},
) where {T}
    linear_solver = get(dest.options, :linear_solver, DirectSolver)
    if T != scsint_t(linear_solver)
        cache = MOI.Utilities.UniversalFallback(
            OptimizerCache{scsint_t(linear_solver)}(),
        )
        MOI.copy_to(cache, src)
        return MOI.optimize!(dest, cache)
    end
    # The real stuff starts here.
    MOI.empty!(dest)
    index_map = MOI.Utilities.identity_index_map(src)
    Ab = src.model.constraints
    A = Ab.coefficients

    for (F, S) in keys(src.constraints)
        throw(MOI.UnsupportedConstraint{F,S}())
    end

    model_attributes = MOI.get(src, MOI.ListOfModelAttributesSet())
    max_sense = false
    obj_attr = nothing
    for attr in model_attributes
        if attr == MOI.ObjectiveSense()
            max_sense = MOI.get(src, attr) == MOI.MAX_SENSE
        elseif attr == MOI.Name()
            continue  # This can be skipped without consequence
        elseif attr isa MOI.ObjectiveFunction
            obj_attr = attr
        else
            throw(MOI.UnsupportedAttribute(attr))
        end
    end
    objective_constant, c, P = 0.0, zeros(A.n), SparseArrays.spzeros(A.n, A.n)
    if obj_attr == MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Cdouble}}()
        obj = MOI.get(src, obj_attr)
        objective_constant = MOI.constant(obj)
        for term in obj.terms
            c[term.variable.value] += (max_sense ? -1 : 1) * term.coefficient
        end
    elseif obj_attr ==
           MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{Cdouble}}()
        obj = MOI.get(src, obj_attr)
        objective_constant = MOI.constant(obj)
        scale = max_sense ? -1 : 1
        for term in obj.affine_terms
            c[term.variable.value] += scale * term.coefficient
        end
        nnz = length(obj.quadratic_terms)
        I, J, V = zeros(Int, nnz), zeros(Int, nnz), zeros(Cdouble, nnz)
        for (i, qterm) in enumerate(obj.quadratic_terms)
            I[i] = min(qterm.variable_1.value, qterm.variable_2.value)
            J[i] = max(qterm.variable_1.value, qterm.variable_2.value)
            V[i] = scale * qterm.coefficient
        end
        P = SparseArrays.sparse(I, J, V, A.n, A.n)
    elseif obj_attr !== nothing
        throw(MOI.UnsupportedAttribute(obj_attr))
    end

    # `model.primal` contains the result of the previous optimization.
    # It is used as a warm start if its length is the same, e.g.
    # probably because no variable and/or constraint has been added.
    if A.n != length(dest.sol.primal)
        resize!(dest.sol.primal, A.n)
        fill!(dest.sol.primal, 0.0)
    end
    if A.m != length(dest.sol.dual)
        resize!(dest.sol.dual, A.m)
        fill!(dest.sol.dual, 0.0)
    end
    if A.m != length(dest.sol.slack)
        resize!(dest.sol.slack, A.m)
        fill!(dest.sol.slack, 0.0)
    end
    # Set starting values and throw error for other variable attributes
    has_warm_start = false
    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for attr in MOI.get(src, MOI.ListOfVariableAttributesSet())
        if attr == MOI.VariableName()
            # Skip
        elseif attr == MOI.VariablePrimalStart()
            has_warm_start = true
            for (i, x) in enumerate(vis_src)
                dest.sol.primal[i] = something(MOI.get(src, attr, x), 0.0)
            end
        else
            throw(MOI.UnsupportedAttribute(attr))
        end
    end
    # Set starting values and throw error for other constraint attributes
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        for attr in MOI.get(src, MOI.ListOfConstraintAttributesSet{F,S}())
            if attr == MOI.ConstraintName()
                # Skip
            elseif attr == MOI.ConstraintPrimalStart()
                has_warm_start = true
                for ci in cis_src
                    start = MOI.get(src, attr, ci)
                    if start !== nothing
                        rows = MOI.Utilities.rows(Ab, ci)
                        dest.sol.slack[rows] .= start
                    end
                end
            elseif attr == MOI.ConstraintDualStart()
                has_warm_start = true
                for ci in cis_src
                    start = MOI.get(src, attr, ci)
                    if start !== nothing
                        rows = MOI.Utilities.rows(Ab, ci)
                        dest.sol.dual[rows] .= start
                    end
                end
            else
                throw(MOI.UnsupportedAttribute(attr))
            end
        end
    end
    options = copy(dest.options)
    if dest.silent
        options[:verbose] = 0
    end
    dest.cones = deepcopy(Ab.sets)
    sol = scs_solve(
        linear_solver,
        A.m,
        A.n,
        A,
        P,
        Ab.constants.b,
        c,
        Ab.sets.num_rows[1],
        Ab.sets.num_rows[2] - Ab.sets.num_rows[1],
        Float64[], # # placeholder: bu
        Float64[], # # placeholder: bl
        _map_sets(MOI.dimension, T, Ab, MOI.SecondOrderCone),
        _map_sets(MOI.side_dimension, T, Ab, ScaledPSDCone),
        div(Ab.sets.num_rows[5] - Ab.sets.num_rows[4], 3),
        div(Ab.sets.num_rows[6] - Ab.sets.num_rows[5], 3),
        vcat(
            _map_sets(set -> set.exponent, Float64, Ab, MOI.PowerCone{Cdouble}),
            _map_sets(
                set -> -set.exponent,
                Float64,
                Ab,
                MOI.DualPowerCone{Cdouble},
            ),
        ),
        dest.sol.primal,
        dest.sol.dual,
        dest.sol.slack;
        warm_start = has_warm_start,
        options...,
    )
    # If the solution is an infeasibility certificate, the objective values are
    # left as infinite, not the value corresponding to the ray.
    if !isfinite(sol.info.dobj)
        sol.info.dobj = -Ab.constants.b' * sol.y
    end
    if !isfinite(sol.info.pobj)
        sol.info.pobj = c' * sol.x
    end
    if sol.ret_val == -4  # SCS_FAILED
        dest.sol = MOISolution(
            sol.x,
            sol.y,
            sol.s,
            sol.ret_val,
            "",
            (max_sense ? -1 : 1) * sol.info.pobj,
            (max_sense ? -1 : 1) * sol.info.dobj,
            objective_constant,
            NaN,
            -1,
        )
    else
        dest.sol = MOISolution(
            sol.x,
            sol.y,
            sol.s,
            sol.ret_val,
            raw_status(sol.info),
            (max_sense ? -1 : 1) * sol.info.pobj,
            (max_sense ? -1 : 1) * sol.info.dobj,
            objective_constant,
            (sol.info.setup_time + sol.info.solve_time) / 1000,
            sol.info.iter,
        )
    end
    return index_map, false
end

function MOI.optimize!(dest::Optimizer, src::MOI.ModelLike)
    linear_solver = get(dest.options, :linear_solver, DirectSolver)
    T = scsint_t(linear_solver)
    cache = MOI.Utilities.UniversalFallback(OptimizerCache{T}())
    index_map = MOI.copy_to(cache, src)
    MOI.optimize!(dest, cache)
    return index_map, false
end

function MOI.get(optimizer::Optimizer, ::MOI.SolveTimeSec)
    return optimizer.sol.solve_time_sec
end

function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
    return optimizer.sol.raw_status
end

"""
    ADMMIterations()

The number of ADMM iterations completed during the solve.
"""
struct ADMMIterations <: MOI.AbstractModelAttribute end

MOI.is_set_by_optimize(::ADMMIterations) = true

function MOI.get(optimizer::Optimizer, ::ADMMIterations)
    return optimizer.sol.iterations
end

# Implements getter for result value and statuses
# SCS returns one of the following integers:
# -7 SCS_INFEASIBLE_INACCURATE
# -6 SCS_UNBOUNDED_INACCURATE
# -5 SCS_SIGINT
# -4 SCS_FAILED
# -3 SCS_INDETERMINATE
# -2 SCS_INFEASIBLE  : primal infeasible, dual unbounded
# -1 SCS_UNBOUNDED   : primal unbounded, dual infeasible
#  0 SCS_UNFINISHED  : never returned, used as placeholder
#  1 SCS_SOLVED
#  2 SCS_SOLVED_INACCURATE
function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    s = optimizer.sol.ret_val
    @assert -7 <= s <= 2
    if s == -7
        return MOI.ALMOST_INFEASIBLE
    elseif s == -6
        return MOI.ALMOST_DUAL_INFEASIBLE
    elseif s == 2
        if occursin("reached time_limit_secs", optimizer.sol.raw_status)
            return MOI.TIME_LIMIT
        elseif occursin("reached max_iters", optimizer.sol.raw_status)
            return MOI.ITERATION_LIMIT
        else
            return MOI.ALMOST_OPTIMAL
        end
    elseif s == -5
        return MOI.INTERRUPTED
    elseif s == -4
        return MOI.INVALID_MODEL
    elseif s == -3
        return MOI.SLOW_PROGRESS
    elseif s == -2
        return MOI.INFEASIBLE
    elseif s == -1
        return MOI.DUAL_INFEASIBLE
    elseif s == 1
        return MOI.OPTIMAL
    else
        @assert s == 0
        return MOI.OPTIMIZE_NOT_CALLED
    end
end

function MOI.get(optimizer::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    value = optimizer.sol.objective_value
    if !MOI.Utilities.is_ray(MOI.get(optimizer, MOI.PrimalStatus()))
        value += optimizer.sol.objective_constant
    end
    return value
end

function MOI.get(optimizer::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    value = optimizer.sol.dual_objective_value
    if !MOI.Utilities.is_ray(MOI.get(optimizer, MOI.DualStatus()))
        value += optimizer.sol.objective_constant
    end
    return value
end

function MOI.get(optimizer::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    elseif optimizer.sol.ret_val in (-3, 1, 2)
        return MOI.FEASIBLE_POINT
    elseif optimizer.sol.ret_val in (-6, -1)
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif optimizer.sol.ret_val == -4
        return MOI.NO_SOLUTION
    end
    return MOI.INFEASIBLE_POINT
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.VariablePrimal,
    vi::MOI.VariableIndex,
)
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.sol.primal[vi.value]
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{F,S},
) where {F,S}
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.sol.slack[MOI.Utilities.rows(optimizer.cones, ci)]
end

function MOI.get(optimizer::Optimizer, attr::MOI.DualStatus)
    if attr.result_index > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    elseif optimizer.sol.ret_val in (-3, 1, 2)
        return MOI.FEASIBLE_POINT
    elseif optimizer.sol.ret_val in (-7, -2)
        return MOI.INFEASIBILITY_CERTIFICATE
    elseif optimizer.sol.ret_val == -4
        return MOI.NO_SOLUTION
    end
    return MOI.INFEASIBLE_POINT
end

function MOI.get(
    optimizer::Optimizer,
    attr::MOI.ConstraintDual,
    ci::MOI.ConstraintIndex{F,S},
) where {F,S}
    MOI.check_result_index_bounds(optimizer, attr)
    return optimizer.sol.dual[MOI.Utilities.rows(optimizer.cones, ci)]
end

MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = 1
