#############################################################################
# SCS.jl
# Wrapper around the SCS solver https://github.com/cvxgrp/scs
# See https://github.com/karanveerm/SCS.jl/
#############################################################################
# SCSSolverInterface.jl
# MathProgBase.jl interface for the SCS.jl solver wrapper
#############################################################################

require(joinpath(Pkg.dir("MathProgBase"), "src", "MathProgSolverInterface.jl"))
importall MathProgSolverInterface
import Base.convert

function convert(x::Type{Int64}, y::UnitRange{Int64})
    if length(y) == 1
        return y[1]
    else
        Base.convert(x, y)
    end
end
#############################################################################
# Define the MPB Solver and Model objects
export SCSSolver
immutable SCSSolver <: AbstractMathProgSolver
end

type SCSMathProgModel <: AbstractMathProgModel
    m::Int64                            # Number of constraints
    n::Int64                            # Number of variables
    A::SparseMatrixCSC{Float64,Int}     # The A matrix (equalities)
    b::Vector{Float64}                  # RHS
    c::Vector{Float64}                  # The objective coeffs (always min)
    f::Int64                            # number of zero cones
    l::Int64                            # number of linear cones { x | x >= 0}
    q::Array{Int64,}                    # Array of SOC sizes
    qsize::Int64                        # Length of q
    s::Array{Int64,}                    # Array of SDC sizes
    ssize::Int64                        # Length of s
    ep::Int64                           # Number of primal exponential cones
    ed::Int64                           # Number of dual exponential cones
    orig_sense::Symbol                  # Original objective sense
    # Post-solve
    solve_stat::Symbol
    obj_val::Float64
    primal_sol::Vector{Float64}
    dual_sol::Vector{Float64}
    slack::Vector{Float64}
    fwd_map::Vector{Int}                # To reorder solution if we solved
end                                     # using the conic interface

SCSMathProgModel() = SCSMathProgModel(0, 0, spzeros(0, 0), Int[], Int[],
                                      0, 0, Int[], 0, Int[], 0, 0, 0,
                                      :Min, :NotSolved, 0.0, Float64[], Float64[],
                                      Float64[],Int[])

#############################################################################
# Begin implementation of the MPB low-level interface
# Implements
# - model
# - loadproblem!
# - optimize!
# - status
# http://mathprogbasejl.readthedocs.org/en/latest/lowlevel.html

model(s::SCSSolver) = SCSMathProgModel()

# Loads the provided problem data to set up the linear programming problem:
# min c'x
# st  lb <= Ax <= ub
#      l <=  x <= u
# where sense = :Min or :Max
function loadproblem!(m::SCSMathProgModel, A, collb, colub, obj, rowlb, rowub, sense)
    (nvar = length(collb)) == length(colub) || error("Unequal lengths for column bounds")
    (nrow = length(rowlb)) == length(rowub) || error("Unequal lengths for row bounds")

    # Turn variable bounds into constraints
    # Inefficient, because keeps allocating memory!
    # Would need to batch, get tricky...
    for j = 1:nvar
        if collb[j] != -Inf
            # Variable has lower bound
            newrow = zeros(1, nvar)
            newrow[j] = -1.0
            A = vcat(A, newrow)
            rowlb = vcat(rowlb, -Inf)
            rowub = vcat(rowub, -collb[j])
            nrow += 1
        end
        if colub[j] != +Inf
            # Variable has upper bound
            newrow = zeros(1, nvar)
            newrow[j] = 1.0
            A = vcat(A, newrow)
            rowlb = vcat(rowlb, -Inf)
            rowub = vcat(rowub, colub[j])
            nrow += 1
        end
    end

    eqidx   = Int[]      # Equality row indicies
    ineqidx = Int[]      # Inequality row indicies
    eqbnd   = Float64[]  # Bounds for equality rows
    ineqbnd = Float64[]  # Bounds for inequality row
    for it in 1:nrow
        # Equality constraint
        if rowlb[it] == rowub[it]
            push!(eqidx, it)
            push!(eqbnd, rowlb[it])
        # Range constraint - not supported
        elseif rowlb[it] != -Inf && rowub[it] != Inf
            error("Ranged constraints unsupported!")
        # Less-than constraint
        elseif rowlb[it] == -Inf
            push!(ineqidx, it)
            push!(ineqbnd, rowub[it])
        # Greater-than constraint - flip sign so only have <= constraints
        else
            push!(ineqidx, it)
            push!(ineqbnd, -rowlb[it])
            A[it,:] *= -1 # flip signs so we have Ax<=b
        end
    end

    m.n         = nvar                              # Number of variables
    m.m         = length(ineqidx) + length(eqidx)   # Number of inequalities Gx <=_K h
    m.A         = sparse([A[eqidx,:]; A[ineqidx,:]])
    m.b         = [eqbnd; ineqbnd]
    m.c         = (sense == :Max) ? obj * -1 : obj[:]
                                        # The objective coeffs (always min)
    m.f         = length(eqidx)
    m.l         = length(ineqidx)
    m.q         = Int64[]
    m.qsize     = 0
    m.s         = Int64[]
    m.ssize     = 0
    m.ep        = 0
    m.ed        = 0
    m.orig_sense = sense                # Original objective sense
    m.fwd_map   = [1:nvar]              # Identity mapping
end

function optimize!(m::SCSMathProgModel)
  solution = SCS_solve(m.m, m.n, m.A, m.b, m.c, m.f, m.l, m.q, m.qsize,
      m.s, m.ssize, m.ep, m.ed)

  m.solve_stat = solution.status
  m.primal_sol = solution.x[m.fwd_map]
  # TODO: Do we need to do some sort of mapping for the dual solution?
  m.dual_sol = solution.y
  m.slack = solution.s
  m.obj_val = dot(m.c, m.primal_sol) * (m.orig_sense == :Max ? -1 : +1)
end

status(m::SCSMathProgModel) = m.solve_stat
getobjval(m::SCSMathProgModel) = m.obj_val
getsolution(m::SCSMathProgModel) = m.primal_sol

#############################################################################
# Begin implementation of the MPB conic interface
# Implements
# - loadconicproblem!
# http://mathprogbasejl.readthedocs.org/en/latest/conic.html

# TODO: This is the world's worst implementation. REDO tomorrow
function orderconesforscs(A, b, cones)
    # Order the cones as
    # free cones
    # zero cones
    # linear cones
    # SOC cones
    # SDP cones
    # Exp primal cones
    # Exp dual cones
    scs_A = nothing

    num_vars = 0
    for (cone, idxs) in cones
       num_vars += length(idxs)
    end

    fwd_map = Array(Int, num_vars)
    multipliers = ones(num_vars, 1)

    k = 1
    for (cone, idxs) in cones
        if cone == :Free
            if scs_A == nothing
                scs_A = A[:, idxs]
            else
                scs_A = [scs_A A[:, idxs]]
            end
            fwd_map[idxs] = k:k + length(idxs) - 1
            k += length(idxs)
        end
    end

    for (cone, idxs) in cones
        if cone == :Zero
            if scs_A == nothing
                scs_A = A[:, idxs]
            else
                scs_A = [scs_A A[:, idxs]]
            end
            fwd_map[idxs] = k:k + length(idxs) - 1
            k += length(idxs)
        end
    end

    for (cone, idxs) in cones
        if cone == :NonNeg
            if scs_A == nothing
                scs_A = A[:, idxs]
            else
                scs_A = [scs_A A[:, idxs]]
            end
            fwd_map[idxs] = k:k + length(idxs) - 1
            k += length(idxs)
        end
    end
    for (cone, idxs) in cones
        if cone == :NonPos
            if scs_A == nothing
                scs_A = -A[:, idxs]
            else
                scs_A = [scs_A A[:, idxs]]
            end
            fwd_map[idxs] = k:k + length(idxs) - 1
            multipliers[k:k + length(idxs) - 1] = -1
            k += length(idxs)
        end
    end

    for (cone, idxs) in cones
        if cone == :SOC
            if scs_A == nothing
                scs_A = A[:, idxs]
            else
                scs_A = [scs_A A[:, idxs]]
            end
            fwd_map[idxs] = k:k + length(idxs) - 1
            k += length(idxs)
        end
    end

    for (cone, idxs) in cones
        if cone == :SDP
            if scs_A == nothing
                scs_A = A[:, idxs]
            else
                scs_A = [scs_A A[:, idxs]]
            end
            fwd_map[idxs] = k:k + length(idxs) - 1
            k += length(idxs)
        end
    end

    for (cone, idxs) in cones
        if cone == :ExpPrimal
            if scs_A == nothing
                scs_A = A[:, idxs]
            else
                scs_A = [scs_A A[:, idxs]]
            end
            fwd_map[idxs] = k:k + length(idxs) - 1
            k += length(idxs)
        end
    end

    for (cone, idxs) in cones
        if cone == :ExpDual
            if scs_A == nothing
                scs_A = A[:, idxs]
            else
                scs_A = [scs_A A[:, idxs]]
            end
            fwd_map[idxs] = k:k + length(idxs) - 1
            k += length(idxs)
        end
    end

    return scs_A, b[:], cones, fwd_map, multipliers
end

# TODO: This is the world's worst implementation. REDO tomorrow
function loadconicproblem!(model::SCSMathProgModel, c, A, b, cones)
    # TODO (if it matters): make this more efficient for sparse A

    # We don't support SOCRotated
    # TODO: We should support SOCRotated
    bad_cones = [:SOCRotated]
    for cone_vars in cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end

    scs_A, scs_b, cones, fwd_map, multipliers = orderconesforscs(A, b, cones)

    num_free = 0
    for (cone, idxs) in cones
        if cone == :Free
            num_free += length(idxs)
        end
    end
    multipliers = multipliers[num_free + 1 : end]
    m, n = size(scs_A)

    rows_G = n - num_free
    G = [zeros(rows_G, num_free) -diagm(multipliers)]

    scs_A = [scs_A; G]
    scs_b = [b; zeros(rows_G, 1)]

    # num zero cones
    f = m
    for (cone, idxs) in cones
        if cone == :Zero
            f += length(idxs)
        end
    end

    # num linear cones
    l = 0
    for (cone, idxs) in cones
        if cone == :NonNeg
            l += length(idxs)
        end
    end

    for (cone, idxs) in cones
        if cone == :NonPos
            l += length(idxs)
        end
    end

    q = Int64[]
    qsize = 0
    for (cone, idxs) in cones
        if cone == :SOC
            qsize += 1
            push!(q, length(idxs))
        end
    end

    s = Int64[]
    ssize = 0
    for (cone, idxs) in cones
        if cone == :SDP
            # n must be a square integer
            n = length(idxs)
            try
                sqrt_n = convert(Int, sqrt(n))
            catch
                error("number of SDP variables must be square")
            end
            ssize += 1
            push!(s, sqrt_n)
        end
    end

    ep = 0
    for (cone, idxs) in cones
        if cone == :ExpPrimal
            ep += length(idxs) / 3
        end
    end

    ed = 0
    for (cone, idxs) in cones
        if cone == :ExpDual
            ed += length(idxs) / 3
        end
    end

    model.n         = n
    model.m         = m + rows_G
    model.A         = scs_A
    model.b         = scs_b[:]
    model.c         = c[:]
    model.q         = q
    model.qsize     = qsize
    model.s         = s
    model.ssize     = ssize
    model.ep        = ep
    model.ed        = ed
    model.orig_sense = :Min
    model.f         = f
    model.l         = l
    model.fwd_map   = fwd_map
end
