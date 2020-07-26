using MathOptInterface
const MOI = MathOptInterface
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const MOIU = MOI.Utilities

"""
    ADMMIterations()

The number of ADMM iterations completed during the solve.
"""
struct ADMMIterations <: MOI.AbstractModelAttribute
end

MOI.is_set_by_optimize(::ADMMIterations) = true

mutable struct MOISolution
    ret_val::Int
    raw_status::String
    primal::Vector{Float64}
    dual::Vector{Float64}
    slack::Vector{Float64}
    objective_value::Float64
    dual_objective_value::Float64
    objective_constant::Float64
    solve_time_sec::Float64
    iterations::Int
end
MOISolution() = MOISolution(0, # SCS_UNFINISHED
                            "", Float64[], Float64[], Float64[], NaN, NaN, NaN,
                            0.0, 0)

# Used to build the data with allocate-load during `copy_to`.
# When `optimize!` is called, a the data is passed to SCS
# using `SCS_solve` and the `ModelData` struct is discarded
mutable struct ModelData{T<:SCSInt}
    m::Int # Number of rows/constraints
    n::Int # Number of cols/variables
    colptr::Vector{T}
    rowval::Vector{T}
    nzval::Vector{Float64}
    b::Vector{Float64} # constants
    objective_constant::Float64 # The objective is min c'x + objective_constant
    c::Vector{Float64}
    function ModelData{T}(n) where T
        data = new()
        data.n = n
        data.colptr = zeros(T, n + 1)
        data.objective_constant = 0.0
        data.c = zeros(n)
        return data
    end
end
function _allocate_nzvals(data::ModelData{T}) where T
    for i in 3:length(data.colptr)
        data.colptr[i] += data.colptr[i - 1]
    end
    data.rowval = Vector{T}(undef, data.colptr[end])
    data.nzval = Vector{Float64}(undef, data.colptr[end])
end
function _managed_matrix(data::ModelData{T}) where T
    for i in length(data.colptr):-1:2
        data.colptr[i] = data.colptr[i - 1]
    end
    data.colptr[1] = 0
    return ManagedSCSMatrix{T}(data.m, data.n, data.nzval, data.rowval, data.colptr)
end

# This is tied to SCS's internal representation
mutable struct ConeData
    f::Int # number of linear equality constraints
    l::Int # length of LP cone
    q::Int # length of SOC cone
    qa::Vector{Int} # array of second-order cone constraints
    s::Int # length of SD cone
    sa::Vector{Int} # array of semi-definite constraints
    ep::Int # number of primal exponential cone triples
    ed::Int # number of dual exponential cone triples
    p::Vector{Float64} # array of power cone params
    nrows::Dict{Int, Int} # The number of rows of each vector sets, this is used by `constrrows` to recover the number of rows used by a constraint when getting `ConstraintPrimal` or `ConstraintDual`
    function ConeData()
        new(0, 0, 0, Int[], 0, Int[], 0, 0, Float64[], Dict{Int, Int}())
    end
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    cone::ConeData
    maxsense::Bool
    data::Union{Nothing, ModelData{Int}, ModelData{Int32}} # only non-Void between MOI.copy_to and MOI.optimize!
    sol::MOISolution
    silent::Bool
    options::Dict{Symbol, Any}
    function Optimizer(; kwargs...)
        optimizer = new(ConeData(), false, nothing, MOISolution(), false,
                        Dict{Symbol, Any}())
        for (key, value) in kwargs
            MOI.set(optimizer, MOI.RawParameter(String(key)), value)
        end
        return optimizer
    end

end

MOI.get(::Optimizer, ::MOI.SolverName) = "SCS"

function MOI.set(optimizer::Optimizer, param::MOI.RawParameter, value)
    # TODO(odow): remove warning in future version.
    if !(param.name isa String)
        Base.depwarn(
            "passing `$(param.name)` to `MOI.RawParameter` as type " *
            "`$(typeof(param.name))` is deprecated. Use a string instead.",
            Symbol("MOI.set")
        )
    end
    optimizer.options[Symbol(param.name)] = value
end
function MOI.get(optimizer::Optimizer, param::MOI.RawParameter)
    # TODO(odow): remove warning in future version.
    if !(param.name isa String)
        Base.depwarn(
            "passing $(param.name) to `MOI.RawParameter` as type " *
            "$(typeof(param.name)) is deprecated. Use a string instead.",
            Symbol("MOI.get")
        )
    end
    return optimizer.options[Symbol(param.name)]
end

MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.set(optimizer::Optimizer, ::MOI.Silent, value::Bool)
    optimizer.silent = value
end
MOI.get(optimizer::Optimizer, ::MOI.Silent) = optimizer.silent

function MOI.is_empty(optimizer::Optimizer)
    !optimizer.maxsense && optimizer.data === nothing
end
function MOI.empty!(optimizer::Optimizer)
    optimizer.maxsense = false
    optimizer.data = nothing # It should already be nothing except if an error is thrown inside copy_to
    optimizer.sol.ret_val = 0
end

function MOI.supports(::Optimizer,
                      ::Union{MOI.ObjectiveSense,
                              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}})
    return true
end

function MOI.supports(::Optimizer, ::MOI.VariablePrimalStart,
                      ::Type{MOI.VariableIndex})
    return true
end

function MOI.supports(::Optimizer,
                      ::Union{MOI.ConstraintPrimalStart,
                              MOI.ConstraintDualStart},
                      ::Type{<:MOI.ConstraintIndex})
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{<:MOI.VectorAffineFunction{Float64}},
    ::Type{<:Union{MOI.Zeros, MOI.Nonnegatives, MOI.SecondOrderCone,
                   MOI.ExponentialCone, MOI.DualExponentialCone,
                   MOI.PositiveSemidefiniteConeTriangle,
                   MOI.PowerCone{Float64}, MOI.DualPowerCone{Float64}}})
    return true
end

const CONE_TYPES = [MOI.Zeros, MOI.Nonnegatives, MOI.SecondOrderCone, MOI.PositiveSemidefiniteConeTriangle, MOI.ExponentialCone, MOI.DualExponentialCone, MOI.PowerCone{Float64}, MOI.DualPowerCone{Float64}]

function _set_index(FS)
    F, S = FS
    i = findfirst(s -> s == S, CONE_TYPES)
    if i === nothing || F != MOI.VectorAffineFunction{Float64}
        throw(MOI.UnsupportedConstraint{F, S}())
    end
    return (i, S)
end

function _allocate_con(optimizer::Optimizer, src, indexmap, ci)
    # TODO use `CanonicalConstraintFunction`
    func = MOI.Utilities.canonical(MOI.get(src, MOI.ConstraintFunction(), ci))
    set = MOI.get(src, MOI.ConstraintSet(), ci)
    ci_dest = _allocate_con(optimizer, indexmap, func, set)
    indexmap[ci] = ci_dest
    return ci_dest, func, set
end

function _allocate_constraints(optimizer::Optimizer, src, indexmap, ::Type{S}) where S
    cis = MOI.get(src, MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{Float64}, S}())
    return map(cis) do ci
        _allocate_con(optimizer, src, indexmap, ci)
    end
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; copy_names::Bool=true)
    MOI.empty!(dest)

    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    idxmap = MOIU.IndexMap()

    constraint_types = MOI.get(src, MOI.ListOfConstraints())
    constraints = map(_set_index, constraint_types)
    sort!(constraints, by = iS -> iS[1])

    _allocate_variables(dest, vis_src, idxmap)

    # Allocate variable attributes
    MOIU.pass_attributes(dest, src, copy_names, idxmap, vis_src, MOIU.allocate)

    # Allocate model attributes
    MOIU.pass_attributes(dest, src, copy_names, idxmap, MOIU.allocate)

    # Allocate constraints
    ci_func_sets = map(constraints) do iS
        _allocate_constraints(dest, src, idxmap, iS[2])
    end

    # Load variables
    MOIU.load_variables(dest, length(vis_src))

    # Load variable attributes
    MOIU.pass_attributes(dest, src, copy_names, idxmap, vis_src, MOIU.load)

    # Load model attributes
    MOIU.pass_attributes(dest, src, copy_names, idxmap, MOIU.load)

    # Load constraints
    for ci_func_set in ci_func_sets
        _load_constraints(dest, idxmap, ci_func_set)
    end

    return idxmap

end

using SparseArrays

# Computes cone dimensions
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractFunction, MOI.Zeros})
    return ci.value
end
#_allocate_constraint: Allocate indices for the constraint `f`-in-`s`
# using information in `cone` and then update `cone`
function _allocate_constraint(cone::ConeData, f, s::MOI.Zeros)
    ci = cone.f
    cone.f += MOI.dimension(s)
    return ci
end
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractFunction, MOI.Nonnegatives})
    return cone.f + ci.value
end
function _allocate_constraint(cone::ConeData, f, s::MOI.Nonnegatives)
    ci = cone.l
    cone.l += MOI.dimension(s)
    return ci
end
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractFunction, MOI.SecondOrderCone})
    return cone.f + cone.l + ci.value
end
function _allocate_constraint(cone::ConeData, f, s::MOI.SecondOrderCone)
    push!(cone.qa, s.dimension)
    ci = cone.q
    cone.q += MOI.dimension(s)
    return ci
end
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractFunction,
                             MOI.PositiveSemidefiniteConeTriangle})
    return cone.f + cone.l + cone.q + ci.value
end
function _allocate_constraint(cone::ConeData, f,
                              s::MOI.PositiveSemidefiniteConeTriangle)
    push!(cone.sa, s.side_dimension)
    ci = cone.s
    cone.s += MOI.dimension(s)
    return ci
end
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractFunction, MOI.ExponentialCone})
    return cone.f + cone.l + cone.q + cone.s + ci.value
end
function _allocate_constraint(cone::ConeData, f, s::MOI.ExponentialCone)
    ci = 3cone.ep
    cone.ep += 1
    return ci
end
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractFunction, MOI.DualExponentialCone})
    return cone.f + cone.l + cone.q + cone.s + 3cone.ep + ci.value
end
function _allocate_constraint(cone::ConeData, f, s::MOI.DualExponentialCone)
    ci = 3cone.ed
    cone.ed += 1
    return ci
end
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractFunction, MOI.PowerCone{Float64}})
    return cone.f + cone.l + cone.q + cone.s + 3cone.ep + 3cone.ed + ci.value
end
function _allocate_constraint(cone::ConeData, f, s::MOI.PowerCone{Float64})
    ci = length(cone.p)
    push!(cone.p, s.exponent)
    return ci
end
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractFunction, MOI.DualPowerCone{Float64}})
    return cone.f + cone.l + cone.q + cone.s + 3cone.ep + 3cone.ed + ci.value
end
function _allocate_constraint(cone::ConeData, f, s::MOI.DualPowerCone{Float64})
    ci = length(cone.p)
    # SCS' convention: dual cones have a negative exponent.
    push!(cone.p, -s.exponent)
    return ci
end
function constroffset(optimizer::Optimizer, ci::CI)
    return constroffset(optimizer.cone, ci::CI)
end
function allocate_terms(colptr, indexmap, terms)
    for term in terms
        colptr[indexmap[term.scalar_term.variable_index].value + 1] += 1
    end
end
function _allocate_con(optimizer::Optimizer, indexmap, func::MOI.VectorAffineFunction{Float64}, s::S) where S <: MOI.AbstractSet
    allocate_terms(optimizer.data.colptr, indexmap, func.terms)
    return CI{MOI.VectorAffineFunction{Float64}, S}(_allocate_constraint(optimizer.cone, func, s))
end

# Vectorized length for matrix dimension n
sympackedlen(n) = div(n*(n+1), 2)
# Matrix dimension for vectorized length n
sympackeddim(n) = div(isqrt(1+8n) - 1, 2)
trimap(i::Integer, j::Integer) = i < j ? trimap(j, i) : div((i-1)*i, 2) + j
trimapL(i::Integer, j::Integer, n::Integer) = i < j ? trimapL(j, i, n) : i + div((2n-j) * (j-1), 2)
function _sympackedto(x, n, mapfrom, mapto)
    @assert length(x) == sympackedlen(n)
    y = similar(x)
    for i in 1:n, j in 1:i
        y[mapto(i, j)] = x[mapfrom(i, j)]
    end
    y
end
sympackedLtoU(x, n=sympackeddim(length(x))) = _sympackedto(x, n, (i, j) -> trimapL(i, j, n), trimap)
sympackedUtoL(x, n=sympackeddim(length(x))) = _sympackedto(x, n, trimap, (i, j) -> trimapL(i, j, n))

function sympackedUtoLidx(x::AbstractVector{<:Integer}, n)
    y = similar(x)
    map =
    for i in eachindex(y)
        y[i] = map[x[i]]
    end
    y
end


# Scale coefficients depending on rows index
# rows: List of row indices
# coef: List of corresponding coefficients
# d: dimension of set
# rev: if true, we unscale instead (e.g. divide by √2 instead of multiply for PSD cone)
function _scalecoef(rows::AbstractVector{<: Integer}, coef::Vector{Float64}, d::Integer, rev::Bool)
    scaling = rev ? 1 / √2 : 1 * √2
    output = copy(coef)
    for i in 1:length(output)
        if !MOI.Utilities.is_diagonal_vectorized_index(rows[i])
            output[i] *= scaling
        end
    end
    return output
end

# Unscale the coefficients in `coef` with respective rows in `rows` for a set `s`
scalecoef(rows, coef, s) = _scalecoef(rows, coef, MOI.dimension(s), false)
# Unscale the coefficients in `coef` with respective rows in `rows` for a set of type `S` with dimension `d`
unscalecoef(rows, coef, d) = _scalecoef(rows, coef, sympackeddim(d), true)

output_index(t::MOI.VectorAffineTerm) = t.output_index
variable_index_value(t::MOI.ScalarAffineTerm) = t.variable_index.value
variable_index_value(t::MOI.VectorAffineTerm) = variable_index_value(t.scalar_term)
coefficient(t::MOI.ScalarAffineTerm) = t.coefficient
coefficient(t::MOI.VectorAffineTerm) = coefficient(t.scalar_term)
# constrrows: Recover the number of rows used by each constraint.
# When, the set is available, simply use MOI.dimension
constrrows(s::MOI.AbstractVectorSet) = 1:MOI.dimension(s)
# When only the index is available, use the `optimizer.ncone.nrows` field
constrrows(optimizer::Optimizer, ci::CI{<:MOI.AbstractVectorFunction, <:MOI.AbstractVectorSet}) = 1:optimizer.cone.nrows[constroffset(optimizer, ci)]

orderval(val, s) = val
function orderval(val, s::MOI.PositiveSemidefiniteConeTriangle)
    sympackedUtoL(val, s.side_dimension)
end
orderidx(idx, s) = idx
function orderidx(idx, s::MOI.PositiveSemidefiniteConeTriangle)
    sympackedUtoLidx(idx, s.side_dimension)
end

# The SCS format is b - Ax ∈ cone
function _load_terms(colptr, rowval, nzval, indexmap, terms, offset)
    @show terms
    @show offset
    for term in terms
        ptr = colptr[indexmap[term.scalar_term.variable_index].value] += 1
        @show ptr
        rowval[ptr] = offset + term.output_index - 1
        @show rowval[ptr]
        nzval[ptr] = -term.scalar_term.coefficient
        @show nzval[ptr]
    end
end
function _load_constant(b, offset, rows, constant::Function)
    for row in rows
        b[offset + row] = constant(row)
    end
end

function _scale(i, coef)
    if MOI.Utilities.is_diagonal_vectorized_index(i)
        return coef
    else
        return coef * √2
    end
end
function _load_con(data, indexmap, func, set::MOI.PositiveSemidefiniteConeTriangle, offset)
    n = set.side_dimension
    map = sympackedLtoU(1:sympackedlen(n), n)
    function map_term(t::MOI.VectorAffineTerm)
        return MOI.VectorAffineTerm(
            map[t.output_index],
            MOI.ScalarAffineTerm(
                _scale(t.output_index, t.scalar_term.coefficient),
                t.scalar_term.variable_index
            )
        )
    end
    new_func = MOI.VectorAffineFunction{Float64}(MOI.VectorAffineTerm{Float64}[map_term(t) for t in func.terms], func.constants)
    # The rows have been reordered in `map_term` so we need to re-canonicalize to reorder the rows.
    MOI.Utilities.canonicalize!(new_func)
    _load_terms(data.colptr, data.rowval, data.nzval, indexmap, new_func.terms, offset)
    function constant(row)
        i = map[row]
        return _scale(i, func.constants[i])
    end
    _load_constant(data.b, offset, eachindex(func.constants), constant)
end
function _load_con(data, indexmap, func, set, offset)
    _load_terms(data.colptr, data.rowval, data.nzval, indexmap, func.terms, offset)
    _load_constant(data.b, offset, eachindex(func.constants), i -> func.constants[i])
end
function _load_constraint(optimizer::Optimizer, indexmap, ci::MOI.ConstraintIndex, func::MOI.VectorAffineFunction{Float64}, s::MOI.AbstractVectorSet)
    offset = constroffset(optimizer, ci)
    _load_con(optimizer.data, indexmap, func, s, offset)
    optimizer.cone.nrows[offset] = MOI.output_dimension(func)
end
function _load_constraints(optimizer::Optimizer, indexmap, ci_func_sets)
    for (ci, func, set) in ci_func_sets
        _load_constraint(optimizer, indexmap, ci, func, set)
    end
end

function _allocate_variables(optimizer::Optimizer, vis_src, idxmap)
    optimizer.cone = ConeData()
    linear_solver, _ = sanitize_SCS_options(optimizer.options)
    T = scsint_t(linear_solver)
    optimizer.data = ModelData{T}(length(vis_src))
    for (i, vi) in enumerate(vis_src)
        idxmap[vi] = MOI.VariableIndex(i)
    end
    return
end

function MOIU.load_variables(optimizer::Optimizer, nvars::Integer)
    cone = optimizer.cone
    m = cone.f + cone.l + cone.q + cone.s + 3cone.ep + 3cone.ed + 3length(cone.p)
    optimizer.data.m = m
    optimizer.data.b = zeros(m)
    _allocate_nzvals(optimizer.data)
    # `optimizer.sol` contains the result of the previous optimization.
    # It is used as a warm start if its length is the same, e.g.
    # probably because no variable and/or constraint has been added.
    if length(optimizer.sol.primal) != nvars
        optimizer.sol.primal = zeros(nvars)
    end
    @assert length(optimizer.sol.dual) == length(optimizer.sol.slack)
    if length(optimizer.sol.dual) != m
        optimizer.sol.dual = zeros(m)
        optimizer.sol.slack = zeros(m)
    end
end

function MOIU.allocate(::Optimizer, ::MOI.VariablePrimalStart,
                       ::MOI.VariableIndex, ::Union{Nothing, Float64})
end
function MOIU.allocate(::Optimizer, ::MOI.ConstraintPrimalStart,
                       ::MOI.ConstraintIndex,
                       ::Union{Nothing, AbstractVector{Float64}})
end
function MOIU.allocate(::Optimizer, ::MOI.ConstraintDualStart,
                       ::MOI.ConstraintIndex,
                       ::Union{Nothing, AbstractVector{Float64}})
end
function MOIU.allocate(optimizer::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    optimizer.maxsense = sense == MOI.MAX_SENSE
end
function MOIU.allocate(::Optimizer, ::MOI.ObjectiveFunction,
                       ::MOI.Union{MOI.ScalarAffineFunction{Float64}})
end

function MOIU.load(::Optimizer, ::MOI.VariablePrimalStart,
                   ::MOI.VariableIndex, ::Nothing)
end
function MOIU.load(optimizer::Optimizer, ::MOI.VariablePrimalStart,
                   vi::MOI.VariableIndex, value::Float64)
    optimizer.sol.primal[vi.value] = value
end
function MOIU.load(::Optimizer, ::MOI.ConstraintPrimalStart,
                   ::MOI.ConstraintIndex, ::Nothing)
end
function MOIU.load(optimizer::Optimizer, ::MOI.ConstraintPrimalStart,
                   ci::MOI.ConstraintIndex, value)
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    optimizer.sol.slack[offset .+ rows] .= value
end
function MOIU.load(::Optimizer, ::MOI.ConstraintDualStart,
                   ::MOI.ConstraintIndex, ::Nothing)
end
function MOIU.load(optimizer::Optimizer, ::MOI.ConstraintDualStart,
                   ci::MOI.ConstraintIndex, value)
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    optimizer.sol.dual[offset .+ rows] .= value
end
function MOIU.load(::Optimizer, ::MOI.ObjectiveSense, ::MOI.OptimizationSense)
end
function MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveFunction,
                   f::MOI.ScalarAffineFunction)
    c0 = Vector(sparsevec(variable_index_value.(f.terms), coefficient.(f.terms),
                          optimizer.data.n))
    optimizer.data.objective_constant = f.constant
    optimizer.data.c = optimizer.maxsense ? -c0 : c0
    return nothing
end

function MOI.optimize!(optimizer::Optimizer)
    cone = optimizer.cone
    m = optimizer.data.m
    n = optimizer.data.n
    b = optimizer.data.b
    objective_constant = optimizer.data.objective_constant
    c = optimizer.data.c

    options = optimizer.options
    if optimizer.silent
        options = copy(options)
        options[:verbose] = 0
    end

    linear_solver, options = sanitize_SCS_options(options)

    sol = SCS_solve(linear_solver, m, n, _managed_matrix(optimizer.data), b, c,
                    cone.f, cone.l, cone.qa, cone.sa, cone.ep, cone.ed, cone.p,
                    optimizer.sol.primal, optimizer.sol.dual,
                    optimizer.sol.slack; options...)

    ret_val = sol.ret_val
    primal = sol.x
    dual = sol.y
    slack = sol.s
    objective_value = (optimizer.maxsense ? -1 : 1) * sol.info.pobj
    dual_objective_value = (optimizer.maxsense ? -1 : 1) * sol.info.dobj
    solve_time = (sol.info.setupTime + sol.info.solveTime) / 1000
    optimizer.sol = MOISolution(ret_val, raw_status(sol.info), primal, dual,
                                slack, objective_value, dual_objective_value,
                                objective_constant, solve_time, sol.info.iter)
end

function MOI.get(optimizer::Optimizer, ::MOI.SolveTime)
    return optimizer.sol.solve_time_sec
end
function MOI.get(optimizer::Optimizer, ::MOI.RawStatusString)
    return optimizer.sol.raw_status
end

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
    if !MOIU.is_ray(MOI.get(optimizer, MOI.PrimalStatus()))
        value += optimizer.sol.objective_constant
    end
    return value
end
function MOI.get(optimizer::Optimizer, attr::MOI.DualObjectiveValue)
    MOI.check_result_index_bounds(optimizer, attr)
    value = optimizer.sol.dual_objective_value
    if !MOIU.is_ray(MOI.get(optimizer, MOI.DualStatus()))
        value += optimizer.sol.objective_constant
    end
    return value
end

function MOI.get(optimizer::Optimizer, attr::MOI.PrimalStatus)
    if attr.N > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    s = optimizer.sol.ret_val
    if s in (-3, 1, 2)
        MOI.FEASIBLE_POINT
    elseif s in (-6, -1)
        MOI.INFEASIBILITY_CERTIFICATE
    else
        MOI.INFEASIBLE_POINT
    end
end
function MOI.get(optimizer::Optimizer, attr::MOI.VariablePrimal, vi::VI)
    MOI.check_result_index_bounds(optimizer, attr)
    optimizer.sol.primal[vi.value]
end
function MOI.get(optimizer::Optimizer, attr::MOI.VariablePrimal, vi::Vector{VI})
    return MOI.get.(optimizer, attr, vi)
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintPrimal,
                 ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    MOI.check_result_index_bounds(optimizer, attr)
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    primal = optimizer.sol.slack[offset .+ rows]
    if S == MOI.PositiveSemidefiniteConeTriangle
        primal = sympackedLtoU(primal)
        primal = unscalecoef(rows, primal, length(rows))
    end
    return primal
end

function MOI.get(optimizer::Optimizer, attr::MOI.DualStatus)
    if attr.N > MOI.get(optimizer, MOI.ResultCount())
        return MOI.NO_SOLUTION
    end
    s = optimizer.sol.ret_val
    if s in (-3, 1, 2)
        MOI.FEASIBLE_POINT
    elseif s in (-7, -2)
        MOI.INFEASIBILITY_CERTIFICATE
    else
        MOI.INFEASIBLE_POINT
    end
end
function MOI.get(optimizer::Optimizer, attr::MOI.ConstraintDual,
                 ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    MOI.check_result_index_bounds(optimizer, attr)
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    dual = optimizer.sol.dual[offset .+ rows]
    if S == MOI.PositiveSemidefiniteConeTriangle
        dual = sympackedLtoU(dual)
        dual = unscalecoef(rows, dual, length(rows))
    end
    return dual
end

MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = 1
