using MathOptInterface
const MOI = MathOptInterface

MOI.Utilities.@product_of_sets(
    Cones,
    MOI.Zeros,
    MOI.Nonnegatives,
    MOI.SecondOrderCone,
    MOI.PositiveSemidefiniteConeTriangle,
    MOI.ExponentialCone,
    MOI.DualExponentialCone,
    MOI.PowerCone{T},
    MOI.DualPowerCone{T},
)

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
        Vector{Cdouble},
        Cones{Cdouble},
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

mutable struct Optimizer <: MOI.AbstractOptimizer
    cones::Union{Nothing,Cones{Cdouble}}
    sol::MOISolution
    silent::Bool
    options::Dict{Symbol,Any}
    function Optimizer(; kwargs...)
        optimizer = new(nothing, MOISolution(), false, Dict{Symbol,Any}())
        for (key, value) in kwargs
            MOI.set(optimizer, MOI.RawOptimizerAttribute(String(key)), value)
        end
        return optimizer
    end
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
### MOI.AbstractModelAttribute
###

function MOI.supports(
    ::Optimizer,
    ::Union{
        MOI.ObjectiveSense,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
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
            MOI.PositiveSemidefiniteConeTriangle,
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

function MOI.copy_to_and_optimize!(
    dest::Optimizer,
    src::MOI.Utilities.UniversalFallback{OptimizerCache{T}},
) where {T}
    linear_solver = get(dest.options, :linear_solver, DirectSolver)
    if T != scsint_t(linear_solver)
        cache = MOI.Utilities.UniversalFallback(
            OptimizerCache{scsint_t(linear_solver)}(),
        )
        MOI.copy_to(cache, src)
        return MOI.copy_to_and_optimize!(dest, cache)
    end
    # The real stuff starts here.
    MOI.empty!(dest)
    index_map = MOI.Utilities.identity_index_map(src)
    Ab = src.model.constraints
    A = Ab.coefficients

    model_attributes = MOI.get(src, MOI.ListOfModelAttributesSet())
    objective_function_attr =
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Cdouble}}()
    for attr in model_attributes
        if attr != MOI.ObjectiveSense() && attr != objective_function_attr
            throw(MOI.UnsupportedAttribute(attr))
        end
    end
    max_sense = false
    if MOI.ObjectiveSense() in model_attributes
        max_sense = MOI.get(src, MOI.ObjectiveSense()) == MOI.MAX_SENSE
    end
    objective_constant = 0.0
    c = zeros(A.n)
    if objective_function_attr in model_attributes
        obj = MOI.get(src, objective_function_attr)
        objective_constant = MOI.constant(obj)
        for term in obj.terms
            c[term.variable.value] += (max_sense ? -1 : 1) * term.coefficient
        end
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
    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    MOI.Utilities.pass_attributes(dest, src, index_map, vis_src)
    # Set starting values and throw error for other constraint attributes
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        MOI.Utilities.pass_attributes(dest, src, index_map, cis_src)
    end
    options = copy(dest.options)
    if dest.silent
        options[:verbose] = 0
    end
    dest.cones = deepcopy(Ab.sets)
    sol = SCS_solve(
        linear_solver,
        A.m,
        A.n,
        A,
        Ab.constants,
        c,
        Ab.sets.num_rows[1],
        Ab.sets.num_rows[2] - Ab.sets.num_rows[1],
        _map_sets(MOI.dimension, T, Ab, MOI.SecondOrderCone),
        _map_sets(
            MOI.side_dimension,
            T,
            Ab,
            MOI.PositiveSemidefiniteConeTriangle,
        ),
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
        options...,
    )
    dest.sol = MOISolution(
        sol.x,
        sol.y,
        sol.s,
        sol.ret_val,
        raw_status(sol.info),
        (max_sense ? -1 : 1) * sol.info.pobj,
        (max_sense ? -1 : 1) * sol.info.dobj,
        objective_constant,
        (sol.info.setupTime + sol.info.solveTime) / 1000,
        sol.info.iter,
    )
    return index_map
end

function MOI.copy_to_and_optimize!(dest::Optimizer, src::MOI.ModelLike)
    linear_solver = get(dest.options, :linear_solver, DirectSolver)
    T = scsint_t(linear_solver)
    cache = MOI.Utilities.UniversalFallback(OptimizerCache{T}())
    index_map = MOI.copy_to(cache, src)
    MOI.copy_to_and_optimize!(dest, cache)
    return index_map
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
        return MOI.ALMOST_OPTIMAL
    elseif s == -5
        return MOI.INTERRUPTED
    elseif s == -4
        return MOI.NUMERICAL_ERROR
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
