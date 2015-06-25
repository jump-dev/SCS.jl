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
    A::SparseMatrixCSC{Float64,Int}   # The A matrix (equalities)
    b::Vector{Float64}                # RHS
    c::Vector{Float64}                # The objective coeffs (always min)
    f::Int                            # number of zero cones
    l::Int                            # number of linear cones { x | x >= 0}
    q::Array{Int,}                    # Array of SOC sizes
    qsize::Int                        # Length of q
    s::Array{Int,}                    # Array of SDP sizes
    ssize::Int                        # Length of s
    ep::Int                           # Number of primal exponential cones
    ed::Int                           # Number of dual exponential cones
    orig_sense::Symbol                # Original objective sense
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

function setsense!(m::SCSMathProgModel, sns::Symbol)
    if m.orig_sense != sns
        sns == :Min || sns == :Max || error("Unrecognized sense $sns")
        m.orig_sense = sns
        m.c *= -1
    end
    nothing
end

function optimize!(m::SCSMathProgModel)
    solution = SCS_solve(m.m, m.n, m.A, m.b, m.c, m.f, m.l, m.q, m.qsize,
                         m.s, m.ssize, m.ep, m.ed,
                         m.primal_sol, m.dual_sol, m.slack; m.options...)

    m.solve_stat = solution.status
    m.primal_sol = solution.x

    m.dual_sol = solution.y

    # TODO: Get the right slack variables in the right order
    m.slack = solution.s

    m.obj_val = dot(m.c, m.primal_sol) * (m.orig_sense == :Max ? -1 : +1)
end

status(m::SCSMathProgModel) = m.solve_stat
getobjval(m::SCSMathProgModel) = m.obj_val
getsolution(m::SCSMathProgModel) = copy(m.primal_sol)

function invertsdconesize(p)
    return (sqrt(8*p+1) - 1) / 2
end

function isintegertol(n)
    return abs(n - convert(Int, n)) < 1e-16
end

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
        nidx = length(idxs)
        if cone == :NonNeg
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
            isintegertol(invertsdconesize(n)) || error("number of SDP variables must be n*(n+1)/2")
            sqrt_n = convert(Int, invertsdconesize(n));
            push!(sqrt_sdp_sizes, sqrt_n)
        end
    end
    for (cone, idxs) in v_cones
        if cone == :SDP
            nidx = length(idxs)
            A_t = [A_t -sparse(idxs, 1:nidx, ones(nidx), num_vars, nidx)]
            b = [b; zeros(nidx)]
             # n must be a square integer
            isintegertol(invertsdconesize(nidx)) || error("number of SDP variables must be n*(n+1)/2")
            sqrt_n = convert(Int, invertsdconesize(nidx));
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


function loadconicproblem!(model::SCSMathProgModel, c, A::SparseMatrixCSC, b, constr_cones, var_cones)
    # TODO: We should support SOCRotated
    bad_cones = [:SOCRotated]
    for cone_vars in constr_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
    end
    for cone_vars in var_cones
        cone_vars[1] in bad_cones && error("Cone type $(cone_vars[1]) not supported")
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

    # rescale model so that vector inner product on R^(n(n+1)/2) matches matrix inner product on S^n
    # (rescale by default)
    if !((:rescale, false) in model.options)
        rescaleconicproblem!(model)
    end
    return model
end

function rescaleconicproblem!(model::SCSMathProgModel)
    SDPstartidx = model.f + model.l + sum(model.q)
    for iSDP = 1:model.ssize
        nSDP = model.s[iSDP]
        sSDP = int(nSDP*(nSDP+1)/2)
        SDPrange = SDPstartidx + 1 : SDPstartidx + sSDP
        model.b[SDPrange] = rescaleSDP(model.b[SDPrange], nSDP, sqrt(2))
        for i = 1:model.n
            model.A[SDPrange,i] = rescaleSDP(model.A[SDPrange,i], nSDP, sqrt(2))
        end
        SDPstartidx += sSDP
    end
    return model
end

function rescaleSDP(x, n, val)
    # scale off-diagonals by val
    M = zeros(n,n)
    idx = find(tril(ones(n,n)))
    M[idx] = x
    dd = diag(M)
    M = (M - diagm(dd))*val + diagm(dd)
    return M[idx]
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
    # reverse the rescaling of the SDP variables
    if !((:rescale, false) in m.options)
        rescaleconicdual!(m, dual)
    end
    return dual
end

function rescaleconicdual!(model::SCSMathProgModel, dual)
    SDPstartidx = model.f + model.l + sum(model.q)
    for iSDP = 1:model.ssize
        nSDP = model.s[iSDP]
        sSDP = int(nSDP*(nSDP+1)/2)
        SDPrange = SDPstartidx + 1 : SDPstartidx + sSDP
        dual[SDPrange] = rescaleSDP(dual[SDPrange], nSDP, 1/sqrt(2))
    end
    return dual
end

# warmstart
# kwargs can be `primal_sol`, `dual_sol`, and `slack`
function setwarmstart!(m::SCSMathProgModel, primal_sol; kwargs...)
    push!(m.options, (:warm_start, true))
    m.primal_sol = primal_sol
    for (k,v) in kwargs
        setfield!(m, k, v)
    end

    # check sizes to prevent segfaults
    nconstr, nvar = size(m.A)
    length(m.primal_sol) == nvar || (m.primal_sol = zeros(nvar))
    length(m.dual_sol) == nconstr || (m.dual_sol = zeros(nconstr))
    length(m.slack) == nconstr || (m.slack = zeros(nconstr))
    m
end
