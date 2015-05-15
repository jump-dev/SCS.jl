#############################################################################
# SCS.jl
# Wrapper around the SCS solver https://github.com/cvxgrp/scs
# See https://github.com/karanveerm/SCS.jl/
#############################################################################
# SCSSolverInterface.jl
# MathProgBase.jl interface for the SCS.jl solver wrapper
#############################################################################

importall MathProgBase.MathProgSolverInterface
import Base.convert

function convert(x::Type{Int}, y::UnitRange{Int})
    if length(y) == 1
        return y[1]
    else
        error("convert` has no method matching convert(::Type{Int}, ::UnitRange{Int})")
    end
end
#############################################################################
# Define the MPB Solver and Model objects
export SCSSolver
immutable SCSSolver <: AbstractMathProgSolver
    options
end
SCSSolver(;kwargs...) = SCSSolver(kwargs)

type SCSMathProgModel <: AbstractMathProgModel
    m::Int                            # Number of constraints
    n::Int                            # Number of variables
    A::SparseMatrixCSC{Float64,Int}     # The A matrix (equalities)
    b::Vector{Float64}                  # RHS
    c::Vector{Float64}                  # The objective coeffs (always min)
    f::Int                            # number of zero cones
    l::Int                            # number of linear cones { x | x >= 0}
    q::Array{Int,}                    # Array of SOC sizes
    qsize::Int                        # Length of q
    s::Array{Int,}                    # Array of SDC sizes
    ssize::Int                        # Length of s
    ep::Int                           # Number of primal exponential cones
    ed::Int                           # Number of dual exponential cones
    orig_sense::Symbol                  # Original objective sense
    # Post-solve
    solve_stat::Symbol
    obj_val::Float64
    primal_sol::Vector{Float64}
    dual_sol::Vector{Float64}
    slack::Vector{Float64}
    row_map_ind::Vector{Int}
    row_map_type::Vector{Symbol}
    options
end

SCSMathProgModel(;kwargs...) = SCSMathProgModel(0, 0, spzeros(0, 0), Int[], Int[],
                                      0, 0, Int[], 0, Int[], 0, 0, 0,
                                      :Min, :NotSolved, 0.0, Float64[], Float64[],
                                      Float64[], Int[], Symbol[], kwargs)

#############################################################################
# Begin implementation of the MPB low-level interface
# Implements
# - model
# - loadproblem!
# - optimize!
# - status
# http://mathprogbasejl.readthedocs.org/en/latest/lowlevel.html

model(s::SCSSolver) = SCSMathProgModel(;s.options...)

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
    m.q         = Int[]
    m.qsize     = 0
    m.s         = Int[]
    m.ssize     = 0
    m.ep        = 0
    m.ed        = 0
    m.orig_sense = sense                # Original objective sense
end

function optimize!(m::SCSMathProgModel)
    solution = SCS_solve(m.m, m.n, m.A, m.b, m.c, m.f, m.l, m.q, m.qsize,
                         m.s, m.ssize, m.ep, m.ed; m.options...)

    m.solve_stat = solution.status
    m.primal_sol = solution.x

    m.dual_sol = solution.y

    # TODO: Get the right slack variables in the right order
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
# - supportedcones
# http://mathprogbasejl.readthedocs.org/en/latest/conic.html

function orderconesforscs(A_in, b_in, c_cones, v_cones)
    # Order the cones as:
    # Free, Zero, NonNeg (NonPos are converted), SOC, SDP, ExpPrimal, ExpDual
    #
    # Returns:
    # - scs_A (A ordered as needed), b
    # - num_free, num_zero, num_linear, soc_sizes, soc_sizes, sqrt_sdp_sizes,
    # - sqrt_sdp_size, num_expprimal, num_expdual

    m, n = size(A_in)
    A_in_t = A_in'
    A_t = spzeros(n,0)
    b = zeros(0)
    row_map_ind = zeros(Int, length(b_in))
    row_map_type = Array(Symbol, length(b_in))

    # First, count the total number of variables
    num_vars = 0
    for (_, idxs) in v_cones
        num_vars += length(idxs)
    end

    num_free = 0
    zeroidx = Int[]
    nonnegidx = Int[]
    nonposidx = Int[]
    socidx = Int[]
    soc_sizes = Int[]
    new_c_cones = Any[]

    for (cone, idxs) in c_cones
        if cone == :Free
            error("Why are you passing in a free constraint?")
        end
        # merge some cones for efficiency
        if cone == :Zero
            append!(zeroidx, idxs)
        elseif cone == :NonNeg
            append!(nonnegidx, idxs)
        elseif cone == :NonPos
            append!(nonposidx, idxs)
        elseif cone == :SOC
            append!(socidx, idxs)
            push!(soc_sizes, length(idxs))
        else
            push!(new_c_cones, (cone,idxs))
        end
    end
    length(zeroidx) > 0 && push!(new_c_cones, (:Zero, zeroidx))
    length(nonnegidx) > 0 && push!(new_c_cones, (:NonNeg, nonnegidx))
    length(nonposidx) > 0 && push!(new_c_cones, (:NonPos, nonposidx))
    length(socidx) > 0 && push!(new_c_cones, (:SOC, socidx))

    for (cone, idxs) in v_cones
        if cone == :Free
            num_free += length(idxs)
        end
    end

    num_zero = 0
    for (cone, idxs) in new_c_cones
        if cone == :Zero
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]
            num_zero += length(idxs)
        end
    end
    for (cone, idxs) in v_cones
        if cone == :Zero
            nidx = length(idxs)
            A_t = [A_t sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]
            num_zero += nidx
        end
    end

    num_lin = 0
    for (cone, idxs) in new_c_cones
        if cone == :NonNeg
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]
            num_lin += length(idxs)
        elseif cone == :NonPos
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t -A_in_t[:,idxs]]
            b = [b; -b_in[idxs,:]]
            num_lin += length(idxs)
        end
    end
    for (cone, idxs) in v_cones
        if cone == :NonNeg
            nidx = length(idxs)
            A_t = [A_t -sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]
            num_lin += nidx
        elseif cone == :NonPos
            A_t = [A_t sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]
            num_lin += nidx
        end
    end

    for (cone, idxs) in new_c_cones
        if cone == :SOC
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]
        end
    end
    for (cone, idxs) in v_cones
        if cone == :SOC
            nidx = length(idxs)
            A_t = [A_t -sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]
            push!(soc_sizes, nidx)
        end
    end

    sqrt_sdp_sizes = Int[]
    for (cone, idxs) in new_c_cones
        if cone == :SDP
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]
            # n must be a square integer
            n = length(idxs)
            isinteger(sqrt(n)) || error("number of SDP variables must be square")
            sqrt_n = convert(Int, sqrt(n));
            push!(sqrt_sdp_sizes, sqrt_n)
        end
    end
    for (cone, idxs) in v_cones
        if cone == :SDP
            nidx = length(idxs)
            A_t = [A_t -sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]
             # n must be a square integer
            isinteger(sqrt(nidx)) || error("number of SDP variables must be square")
            sqrt_n = convert(Int, sqrt(nidx));
            push!(sqrt_sdp_sizes, sqrt_n)
        end
    end

    num_expprimal = 0
    for (cone, idxs) in new_c_cones
        if cone == :ExpPrimal
            length(idxs) % 3 == 0 ||
                error("Number of ExpPrimal variables must be a multiple of 3")
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]

            num_expprimal += div(length(idxs), 3)
        end
    end
    for (cone, idxs) in v_cones
        if cone == :ExpPrimal
            length(idxs) % 3 == 0 ||
                error("Number of ExpPrimal variables must be a multiple of 3")
            nidx = length(idxs)
            A_t = [A_t -sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]

            num_expprimal += div(length(idxs), 3)
        end
    end

    num_expdual = 0
    for (cone, idxs) in new_c_cones
        if cone == :ExpDual
            row_map_ind[idxs] = size(A_t, 2)+1:size(A_t, 2)+length(idxs)
            row_map_type[idxs] = [cone for i in 1:length(idxs)]

            length(idxs) % 3 == 0 ||
                error("Number of ExpDual variables must be a multiple of 3")
            A_t = [A_t A_in_t[:,idxs]]
            b = [b; b_in[idxs,:]]

            num_expdual += div(length(idxs), 3)
        end
    end
    for (cone, idxs) in v_cones
        if cone == :ExpDual
            length(idxs) % 3 == 0 ||
                error("Number of ExpDual variables must be a multiple of 3")
            nidx = length(idxs)
            A_t = [A_t -sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]

            num_expdual += div(length(idxs), 3)
        end
    end

    return A_t', b, num_free, num_zero, num_lin, soc_sizes,
           sqrt_sdp_sizes, num_expprimal, num_expdual, row_map_ind, row_map_type
end


loadconicproblem!(model::SCSMathProgModel, c, A, b, constr_cones, var_cones) =
    loadconicproblem!(model, c, sparse(A), b, constr_cones, var_cones)


# TODO: At least for Convex.jl, we are guaranteed what A and b look like
# This should allow us to make it far more efficient

# Given an n x n matrix vectorized in the form Ax - b where x is an unknown variable
# Let y = b - Ax and Y = reshape(y, n, n)
# This function then returns a matrix G and vector h
# s.t Gx + s == h, s in the zero cone enforces Y[i, j] == Y[j, i]
function addsymmetryconstraints(A, b)
    size_b = length(b)
    isinteger(sqrt(size_b)) || error("number of SDP variables must be square")
    n = convert(Int, sqrt(size_b));

    # Enforcing Y[i, j] == Y[j, i] is equivalent to enforcing
    # y[idx1] == y[idx2]; idx1 = n * (j - 1) + i; idx2 = n * (i - 1) + j
    # this is equivalent to (A[idx1, :] - A[idx2, :])x == b[idx1] - b[idx2]
    # and we're done

    num_constraints = int(n * (n - 1) / 2)

    G = spzeros(num_constraints, size(A, 2))
    h = zeros(num_constraints, 1)

    idx = 1
    for i = 1:n
        for j = 1:n
            if i >= j
                continue
            end

            idx1 = n * (j - 1) + i
            idx2 = n * (i - 1) + j
            h[idx, 1] = b[idx1] - b[idx2]
            G[idx, :] = A[idx1, :] - A[idx2, :]
            idx += 1
        end
    end
    return G, h
end


function loadconicproblem!(model::SCSMathProgModel, c, A::SparseMatrixCSC, b, constr_cones, var_cones)
    # TODO: We should support SOCRotated
    bad_cones = [:SOCRotated]
    for cone_vars in constr_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end
    for cone_vars in var_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end

    for cone_vars in constr_cones
        # If Ax + s = b, s in SDP cone
        if cone_vars[1] == :SDP
            idxs = cone_vars[2]
            G, h = addsymmetryconstraints(A[idxs, :], b[idxs])
            cur_size = size(A, 1)
            A = [A; G]
            b = [b; h]
            push!(constr_cones, (:Zero, cur_size + 1:size(A, 1)))
        end
    end

    for cone_vars in var_cones
        # If x in SDP cone
        if cone_vars[1] == :SDP
            idxs = cone_vars[2]
            size_x = length(idxs)

            # Basically enforcing x = Ix + 0 is in SDP cone
            G, h = addsymmetryconstraints(eye(size_x), zeros(size_x, 1))
            cur_size = size(A, 1)
            A = [A; G]
            b = [b; h]
            push!(constr_cones, (:Zero, cur_size + 1:size(A, 1)))
        end
    end

    # Convert idxs to an array
    c_cones = [(cone, [idxs...]) for (cone, idxs) in constr_cones]
    v_cones = [(cone, [idxs...]) for (cone, idxs) in var_cones]

    scs_A, scs_b, num_free, f, l, q, s, ep, ed, row_map_ind, row_map_type =
        orderconesforscs(A, b, c_cones, v_cones)

    m, n = size(scs_A)

    model.n             = n
    model.m             = m # + rows_G
    model.A             = scs_A
    model.b             = scs_b[:]
    model.c             = c[:]
    model.q             = q
    model.qsize         = length(q)
    model.s             = s
    model.ssize         = length(s)
    model.ep            = ep
    model.ed            = ed
    model.orig_sense    = :Min
    model.f             = f
    model.l             = l
    model.row_map_ind   = row_map_ind
    model.row_map_type  = row_map_type
    return model
end

supportedcones(s::SCSSolver) = [:Free, :Zero, :NonNeg, :NonPos, :SOC, :SDP, :ExpPrimal, :ExpDual]

function getconicdual(m::SCSMathProgModel)
    # TODO: Why am I flipping signs? Also, do I need to flip the signs for every cone that isnt NonNeg?
    # Flipping signs seems to make it pass the MPB dual tests
    dual = m.dual_sol[m.row_map_ind]
    for i in 1:length(m.row_map_type)
        if m.row_map_type[i] != :NonNeg
            dual[i] = -dual[i]
        end
    end
    return dual
end
