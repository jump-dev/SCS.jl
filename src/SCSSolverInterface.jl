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
    fwd_map::Vector{Int}                # To reorder solution if we solved
end                                   # using the conic interface

SCSMathProgModel() = SCSMathProgModel(0, 0, spzeros(0, 0), Int[], Int[],
                                      0, 0, Int[], 0, Int[], 0, 0, 0,
                                      :Min, :NotSolved, 0.0, Float64[], Float64[],
                                      Int[])

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
    # TODO: There is probably a better way to do this that I can't think of at
    # 12am
    m.A         = sparse([ A[eqidx,:]; A[ineqidx,:] ])
    m.b         = [ eqbnd; ineqbnd ]
    m.c         = (sense == :Max) ? obj * -1 : obj[:]
                                        # The objective coeffs (always min)
    m.q         = Int64[]
    m.qsize     = 0
    m.s         = Int64[]
    m.ssize     = 0
    m.ep        = 0
    m.ed        = 0
    m.orig_sense = sense                # Original objective sense
    m.fwd_map   = [1:nvar]              # Identity mapping

    m.f         = length(eqidx)
    m.l         = length(ineqidx)
end

function optimize!(m::SCSMathProgModel)
  solution = SCS_solve(m.m, m.n, m.A, m.b, m.c, m.f, m.l, m.q, m.qsize,
      m.s, m.ssize, m.ep, m.ed)

  m.solve_stat = solution.status
  m.primal_sol = solution.x
  m.dual_sol = solution.y
  m.obj_val = dot(m.c, m.primal_sol) * (m.orig_sense == :Max ? -1 : +1)
end

status(m::SCSMathProgModel) = m.solve_stat
getobjval(m::SCSMathProgModel) = m.obj_val
getsolution(m::SCSMathProgModel) = m.primal_sol[m.fwd_map]

#############################################################################
# Begin implementation of the MPB conic interface
# Implements
# - loadconicproblem!
# http://mathprogbasejl.readthedocs.org/en/latest/conic.html

function loadconicproblem!(s::SCSMathProgModel, c, A, b, cones)
    # TODO (if it matters): make this more efficient for sparse A

    # We don't support SOCRotated
    # TODO: We should support SOCRotated
    bad_cones = [:SOCRotated]
    for cone_vars in cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end

    G = -speye(n)
    m, n = size(A)

    scs_A = [A; G]
    scs_b = [b; zeros(n, 1)]

    # MathProgBase form             SCS form
    # min c'x                       min c'x
    # st  A x = b                   st  A x + s = b
    #       x in K                      s in K

    # Expand out the cones info
    # The cones can come in any order, so we need to build a mapping
    # from the variable indices in the input to the internal ordering
    # we will use.

    # In the first past we'll just count up the number of variables
    # of each type.
    num_vars = 0
    for (cone_type, idxs) in cones
        num_vars += length(idxs)
    end
    fwd_map = Array(Int,    num_vars)  # Will be used for SOCs
    rev_map = Array(Int,    num_vars)  # Need to restore sol. vec.
    idxcone = Array(Symbol, num_vars)  # We'll uses this for non-SOCs

    # Now build the mapping
    pos = 1
    for (cone, idxs) in cones
        for i in idxs
            fwd_map[i]   = pos   # fwd_map = orig idx -> internal idx
            rev_map[pos] = i     # rev_map = internal idx -> orig idx
            idxcone[pos] = cone
            pos += 1
        end
    end

    # Rearrange data into the internal ordering
    ecos_c = c[rev_map]
    ecos_A = A[:,rev_map]
    ecos_b = b[:]

    # For all variables in the :Zero cone, fix at 0 with an
    # equality constraint. TODO: Don't even include them

    for j = 1:num_vars
        idxcone[j] != :Zero && continue

        new_row    = zeros(1,num_vars)
        new_row[j] = 1.0
        ecos_A     = vcat(ecos_A, new_row)
        ecos_b     = vcat(ecos_b, 0.0)
    end

    # Build G matrix
    # There will be one row for every :NonNeg and :NonPos cone
    # and an additional row for every variable in a :SOC cone
    # Or in other words, everything that isn't a :Free or :Zero
    # gets a row in G and h
    num_G_row = 0
    for j = 1:num_vars
        idxcone[j] == :Free && continue
        idxcone[j] == :Zero && continue
        num_G_row += 1
    end
    ecos_G = zeros(num_G_row,num_vars)
    ecos_h = zeros(num_G_row)

    # First, handle the :NonNeg, :NonPos cases
    num_pos_orth = 0
    G_row = 1
    for j = 1:num_vars
        if idxcone[j] == :NonNeg
            ecos_G[G_row,j] = -1.0
            G_row += 1
            num_pos_orth += 1
        elseif idxcone[j] == :NonPos
            ecos_G[G_row,j] = +1.0
            G_row += 1
            num_pos_orth += 1
        end
    end
    @assert G_row == num_pos_orth + 1
    # Now handle the SOCs
    # The MPB unput form is basically just says a vector of
    # variables (y,x) lives in the SOC  || x || <= y
    # ECOS wants somethings in the form h - Gx in Q so we
    # will prove 0 - Ix \in Q
    num_SOC_cones = 0
    SOC_conedims  = Int[]
    for (cone, idxs) in cones
        cone != :SOC && continue
        # Found a new SOC
        num_SOC_cones += 1
        push!(SOC_conedims, length(idxs))
        # Add the entries (carrying on from pos. orthant)
        for j in idxs
            ecos_G[G_row,fwd_map[j]] = -1.0
            G_row += 1
        end
    end
    @assert G_row == num_G_row + 1

    # Store in the ECOS structure
    m.nvar          = num_vars
    m.nineq         = num_G_row
    m.neq           = length(ecos_b)
    m.npos          = num_pos_orth
    m.ncones        = num_SOC_cones
    m.conedims      = SOC_conedims
    m.G             = ecos_G
    m.A             = ecos_A
    m.c             = ecos_c
    m.orig_sense    = :Min
    m.h             = ecos_h
    m.b             = ecos_b
    m.fwd_map       = fwd_map
end
