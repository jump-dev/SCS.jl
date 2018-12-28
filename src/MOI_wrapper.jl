using MathOptInterface
const MOI = MathOptInterface
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const MOIU = MOI.Utilities

const SF = Union{MOI.SingleVariable, MOI.ScalarAffineFunction{Float64}, MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64}}
const SS = Union{MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives, MOI.SecondOrderCone, MOI.ExponentialCone, MOI.PositiveSemidefiniteConeTriangle, MOI.PowerCone, MOI.DualPowerCone}

mutable struct MOISolution
    ret_val::Int
    primal::Vector{Float64}
    dual::Vector{Float64}
    slack::Vector{Float64}
    objval::Float64
    solve_time::Float64
end
MOISolution() = MOISolution(0, # SCS_UNFINISHED
                            Float64[], Float64[], Float64[], NaN, 0.0)

# Used to build the data with allocate-load during `copy_to`.
# When `optimize!` is called, a the data is passed to SCS
# using `SCS_solve` and the `ModelData` struct is discarded
mutable struct ModelData
    m::Int # Number of rows/constraints
    n::Int # Number of cols/variables
    I::Vector{Int} # List of rows
    J::Vector{Int} # List of cols
    V::Vector{Float64} # List of coefficients
    b::Vector{Float64} # constants
    objconstant::Float64 # The objective is min c'x + objconstant
    c::Vector{Float64}
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
    setconstant::Dict{Int, Float64} # For the constant of EqualTo, LessThan and GreaterThan, they are used for getting the `ConstraintPrimal` as the slack is Ax - b but MOI expects Ax so we need to add the constant b to the slack to get Ax
    nrows::Dict{Int, Int} # The number of rows of each vector sets, this is used by `constrrows` to recover the number of rows used by a constraint when getting `ConstraintPrimal` or `ConstraintDual`
    function ConeData()
        new(0, 0,
            0, Int[],
            0, Int[],
            0, 0, Float64[],
            Dict{Int, Float64}(),
            Dict{Int, Int}())
    end
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    cone::ConeData
    maxsense::Bool
    data::Union{Nothing, ModelData} # only non-Void between MOI.copy_to and MOI.optimize!
    sol::MOISolution
    options
    function Optimizer(; options...)
        new(ConeData(), false, nothing, MOISolution(), options)
    end

end

MOI.get(::Optimizer, ::MOI.SolverName) = "SCS"

function MOI.is_empty(optimizer::Optimizer)
    !optimizer.maxsense && optimizer.data === nothing
end
function MOI.empty!(optimizer::Optimizer)
    optimizer.maxsense = false
    optimizer.data = nothing # It should already be nothing except if an error is thrown inside copy_to
    optimizer.sol.ret_val = 0
end

MOIU.supports_allocate_load(::Optimizer, copy_names::Bool) = !copy_names

function MOI.supports(::Optimizer,
                      ::Union{MOI.ObjectiveSense,
                              MOI.ObjectiveFunction{MOI.SingleVariable},
                              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}},
                              MOI.VariablePrimalStart})
    return true
end

MOI.supports_constraint(::Optimizer, ::Type{<:SF}, ::Type{<:SS}) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    return MOIU.automatic_copy_to(dest, src; kws...)
end

using Compat.SparseArrays

const ZeroCones = Union{MOI.EqualTo, MOI.Zeros}
const LPCones = Union{MOI.GreaterThan, MOI.LessThan, MOI.Nonnegatives, MOI.Nonpositives}

# Computes cone dimensions
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:ZeroCones}) = ci.value
#_allocate_constraint: Allocate indices for the constraint `f`-in-`s` using information in `cone` and then update `cone`
function _allocate_constraint(cone::ConeData, f, s::ZeroCones)
    ci = cone.f
    cone.f += MOI.dimension(s)
    return ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:LPCones}) = cone.f + ci.value
function _allocate_constraint(cone::ConeData, f, s::LPCones)
    ci = cone.l
    cone.l += MOI.dimension(s)
    return ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:MOI.SecondOrderCone}) = cone.f + cone.l + ci.value
function _allocate_constraint(cone::ConeData, f, s::MOI.SecondOrderCone)
    push!(cone.qa, s.dimension)
    ci = cone.q
    cone.q += MOI.dimension(s)
    return ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:MOI.PositiveSemidefiniteConeTriangle}) = cone.f + cone.l + cone.q + ci.value
function _allocate_constraint(cone::ConeData, f, s::MOI.PositiveSemidefiniteConeTriangle)
    push!(cone.sa, s.side_dimension)
    ci = cone.s
    cone.s += MOI.dimension(s)
    return ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:MOI.ExponentialCone}) = cone.f + cone.l + cone.q + cone.s + ci.value
function _allocate_constraint(cone::ConeData, f, s::MOI.ExponentialCone)
    ci = 3cone.ep
    cone.ep += 1
    return ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:MOI.PowerCone}) = cone.f + cone.l + cone.q + cone.s + cone.ep + ci.value
function _allocate_constraint(cone::ConeData, f, s::MOI.PowerCone)
    ci = length(cone.p)
    push!(cone.p, s.exponent)
    return ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, <:MOI.DualPowerCone}) = cone.f + cone.l + cone.q + cone.s + cone.ep + ci.value
function _allocate_constraint(cone::ConeData, f, s::MOI.DualPowerCone)
    ci = length(cone.p)
    # SCS' convention: dual cones have a negative exponent. 
    push!(cone.p, - s.exponent)
    return ci
end
constroffset(optimizer::Optimizer, ci::CI) = constroffset(optimizer.cone, ci::CI)
function MOIU.allocate_constraint(optimizer::Optimizer, f::F, s::S) where {F <: MOI.AbstractFunction, S <: MOI.AbstractSet}
    return CI{F, S}(_allocate_constraint(optimizer.cone, f, s))
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
    map = sympackedLtoU(1:sympackedlen(n), n)
    for i in eachindex(y)
        y[i] = map[x[i]]
    end
    y
end


# Scale coefficients depending on rows index
# rows: List of row indices
# coef: List of corresponding coefficients
# minus: if true, multiply the result by -1
# d: dimension of set
# rev: if true, we unscale instead (e.g. divide by √2 instead of multiply for PSD cone)
_scalecoef(rows, coef, minus, ::Type{<:MOI.AbstractSet}, d, rev) = minus ? -coef : coef
_scalecoef(rows, coef, minus, ::Union{Type{<:MOI.LessThan}, Type{<:MOI.Nonpositives}}, d, rev) = minus ? coef : -coef
function _scalecoef(rows, coef, minus, ::Type{MOI.PositiveSemidefiniteConeTriangle}, d, rev)
    scaling = minus ? -1 : 1
    scaling2 = rev ? scaling / √2 : scaling * √2
    output = copy(coef)
    diagidx = BitSet()
    for i in 1:d
        push!(diagidx, trimap(i, i))
    end
    for i in 1:length(output)
        if rows[i] in diagidx
            output[i] *= scaling
        else
            output[i] *= scaling2
        end
    end
    output
end
# Unscale the coefficients in `coef` with respective rows in `rows` for a set `s` and multiply by `-1` if `minus` is `true`.
scalecoef(rows, coef, minus, s) = _scalecoef(rows, coef, minus, typeof(s), MOI.dimension(s), false)
# Unscale the coefficients in `coef` with respective rows in `rows` for a set of type `S` with dimension `d`
unscalecoef(rows, coef, S::Type{<:MOI.AbstractSet}, d) = _scalecoef(rows, coef, false, S, d, true)
unscalecoef(rows, coef, S::Type{MOI.PositiveSemidefiniteConeTriangle}, d) = _scalecoef(rows, coef, false, S, sympackeddim(d), true)

output_index(t::MOI.VectorAffineTerm) = t.output_index
variable_index_value(t::MOI.ScalarAffineTerm) = t.variable_index.value
variable_index_value(t::MOI.VectorAffineTerm) = variable_index_value(t.scalar_term)
coefficient(t::MOI.ScalarAffineTerm) = t.coefficient
coefficient(t::MOI.VectorAffineTerm) = coefficient(t.scalar_term)
# constrrows: Recover the number of rows used by each constraint.
# When, the set is available, simply use MOI.dimension
constrrows(::MOI.AbstractScalarSet) = 1
constrrows(s::MOI.AbstractVectorSet) = 1:MOI.dimension(s)
# When only the index is available, use the `optimizer.ncone.nrows` field
constrrows(optimizer::Optimizer, ci::CI{<:MOI.AbstractScalarFunction, <:MOI.AbstractScalarSet}) = 1
constrrows(optimizer::Optimizer, ci::CI{<:MOI.AbstractVectorFunction, <:MOI.AbstractVectorSet}) = 1:optimizer.cone.nrows[constroffset(optimizer, ci)]
MOIU.load_constraint(optimizer::Optimizer, ci, f::MOI.SingleVariable, s) = MOIU.load_constraint(optimizer, ci, MOI.ScalarAffineFunction{Float64}(f), s)
function MOIU.load_constraint(optimizer::Optimizer, ci, f::MOI.ScalarAffineFunction, s::MOI.AbstractScalarSet)
    a = sparsevec(variable_index_value.(f.terms), coefficient.(f.terms))
    # sparsevec combines duplicates with + but does not remove zeros created so we call dropzeros!
    dropzeros!(a)
    offset = constroffset(optimizer, ci)
    row = constrrows(s)
    i = offset + row
    # The SCS format is b - Ax ∈ cone
    # so minus=false for b and minus=true for A
    setconstant = MOIU.getconstant(s)
    optimizer.cone.setconstant[offset] = setconstant
    constant = f.constant - setconstant
    optimizer.data.b[i] = scalecoef(row, constant, false, s)
    append!(optimizer.data.I, fill(i, length(a.nzind)))
    append!(optimizer.data.J, a.nzind)
    append!(optimizer.data.V, scalecoef(row, a.nzval, true, s))
end
MOIU.load_constraint(optimizer::Optimizer, ci, f::MOI.VectorOfVariables, s) = MOIU.load_constraint(optimizer, ci, MOI.VectorAffineFunction{Float64}(f), s)
orderval(val, s) = val
function orderval(val, s::MOI.PositiveSemidefiniteConeTriangle)
    sympackedUtoL(val, s.side_dimension)
end
orderidx(idx, s) = idx
function orderidx(idx, s::MOI.PositiveSemidefiniteConeTriangle)
    sympackedUtoLidx(idx, s.side_dimension)
end
function MOIU.load_constraint(optimizer::Optimizer, ci, f::MOI.VectorAffineFunction, s::MOI.AbstractVectorSet)
    A = sparse(output_index.(f.terms), variable_index_value.(f.terms), coefficient.(f.terms))
    # sparse combines duplicates with + but does not remove zeros created so we call dropzeros!
    dropzeros!(A)
    I, J, V = findnz(A)
    offset = constroffset(optimizer, ci)
    rows = constrrows(s)
    optimizer.cone.nrows[offset] = length(rows)
    i = offset .+ rows
    # The SCS format is b - Ax ∈ cone
    # so minus=false for b and minus=true for A
    optimizer.data.b[i] = scalecoef(rows, orderval(f.constants, s), false, s)
    append!(optimizer.data.I, offset .+ orderidx(I, s))
    append!(optimizer.data.J, J)
    append!(optimizer.data.V, scalecoef(I, V, true, s))
end

function MOIU.allocate_variables(optimizer::Optimizer, nvars::Integer)
    optimizer.cone = ConeData()
    VI.(1:nvars)
end

function MOIU.load_variables(optimizer::Optimizer, nvars::Integer)
    cone = optimizer.cone
    m = cone.f + cone.l + cone.q + cone.s + 3cone.ep + cone.ed + 3 * length(cone.p)
    I = Int[]
    J = Int[]
    V = Float64[]
    b = zeros(m)
    c = zeros(nvars)
    optimizer.data = ModelData(m, nvars, I, J, V, b, 0., c)
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
                       ::MOI.VariableIndex, ::Float64)
end
function MOIU.allocate(::Optimizer, ::MOI.ConstraintPrimalStart,
                       ::MOI.ConstraintIndex, ::Float64)
end
function MOIU.allocate(::Optimizer, ::MOI.ConstraintDualStart,
                       ::MOI.ConstraintIndex, ::Float64)
end
function MOIU.allocate(optimizer::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    optimizer.maxsense = sense == MOI.MAX_SENSE
end
function MOIU.allocate(::Optimizer, ::MOI.ObjectiveFunction,
                       ::MOI.Union{MOI.SingleVariable,
                                   MOI.ScalarAffineFunction{Float64}})
end

function MOIU.load(optimizer::Optimizer, ::MOI.VariablePrimalStart,
                   vi::MOI.VariableIndex, value::Float64)
    optimizer.sol.primal[vi.value] = value
end
function MOIU.load(optimizer::Optimizer, ::MOI.ConstraintPrimalStart,
                   ci::MOI.ConstraintIndex, value)
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    optimizer.sol.primal[offset .+ rows] .= value
end
function MOIU.load(optimizer::Optimizer, ::MOI.ConstraintDualStart,
                   ci::MOI.ConstraintIndex, value)
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    optimizer.sol.primal[offset .+ rows] .= value
end
function MOIU.load(::Optimizer, ::MOI.ObjectiveSense, ::MOI.OptimizationSense)
end
function MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveFunction,
                   f::MOI.SingleVariable)
    MOIU.load(optimizer,
              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
              MOI.ScalarAffineFunction{Float64}(f))
end
function MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveFunction,
                   f::MOI.ScalarAffineFunction)
    c0 = Vector(sparsevec(variable_index_value.(f.terms), coefficient.(f.terms),
                          optimizer.data.n))
    optimizer.data.objconstant = f.constant
    optimizer.data.c = optimizer.maxsense ? -c0 : c0
    return nothing
end

function MOI.optimize!(optimizer::Optimizer)
    cone = optimizer.cone
    m = optimizer.data.m
    n = optimizer.data.n
    A = sparse(optimizer.data.I, optimizer.data.J, optimizer.data.V)
    b = optimizer.data.b
    objconstant = optimizer.data.objconstant
    c = optimizer.data.c
    optimizer.data = nothing # Allows GC to free optimizer.data before A is loaded to SCS

    linear_solver, options = sanatize_SCS_options(optimizer.options)

    t = time()
    sol = SCS_solve(linear_solver, m, n, A, b, c,
                    cone.f, cone.l, cone.qa, cone.sa, cone.ep, cone.ed, cone.p,
                    optimizer.sol.primal, optimizer.sol.dual,
                    optimizer.sol.slack; options...)
    solve_time = time() - t

    ret_val = sol.ret_val
    primal = sol.x
    dual = sol.y
    slack = sol.s
    objval = (optimizer.maxsense ? -1 : 1) * dot(c, primal) + objconstant
    optimizer.sol = MOISolution(ret_val, primal, dual, slack, objval, solve_time)
end

function MOI.get(optimizer::Optimizer, ::MOI.SolveTime)
    return optimizer.sol.solve_time
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

MOI.get(optimizer::Optimizer, ::MOI.ObjectiveValue) = optimizer.sol.objval

function MOI.get(optimizer::Optimizer, ::MOI.PrimalStatus)
    s = optimizer.sol.ret_val
    if s in (-3, 1, 2)
        MOI.FEASIBLE_POINT
    elseif s in (-6, -1)
        MOI.INFEASIBILITY_CERTIFICATE
    else
        MOI.INFEASIBLE_POINT
    end
end
function MOI.get(optimizer::Optimizer, ::MOI.VariablePrimal, vi::VI)
    optimizer.sol.primal[vi.value]
end
MOI.get(optimizer::Optimizer, a::MOI.VariablePrimal, vi::Vector{VI}) = MOI.get.(optimizer, a, vi)
_unshift(optimizer::Optimizer, offset, value, s) = value
_unshift(optimizer::Optimizer, offset, value, s::Type{<:MOI.AbstractScalarSet}) = value + optimizer.cone.setconstant[offset]
reorderval(val, s) = val
function reorderval(val, ::Type{MOI.PositiveSemidefiniteConeTriangle})
    sympackedLtoU(val)
end
function MOI.get(optimizer::Optimizer, ::MOI.ConstraintPrimal, ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    _unshift(optimizer, offset, unscalecoef(rows, reorderval(optimizer.sol.slack[offset .+ rows], S), S, length(rows)), S)
end

function MOI.get(optimizer::Optimizer, ::MOI.DualStatus)
    s = optimizer.sol.ret_val
    if s in (-3, 1, 2)
        MOI.FEASIBLE_POINT
    elseif s in (-7, -2)
        MOI.INFEASIBILITY_CERTIFICATE
    else
        MOI.INFEASIBLE_POINT
    end
end
function MOI.get(optimizer::Optimizer, ::MOI.ConstraintDual, ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    unscalecoef(rows, reorderval(optimizer.sol.dual[offset .+ rows], S), S, length(rows))
end

MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = 1
