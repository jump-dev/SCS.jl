#############################################################################
# SCS.jl
# Wrapper around the SCS solver https://github.com/cvxgrp/scs
# See http://github.com/JuliaOpt/SCS.jl
#############################################################################
# test/options.jl
# Tests the ability to pass options
#############################################################################

using Base.Test
using MathProgBase.SolverInterface
using SCS


# The normal test
A = [1.0 1.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 1.0;
     0.0 0.0 1.0 1.0 1.0]
collb = [0.0, 0.0, 0.0, 0.0, 0.0]
obj   = [3.0, 4.0, 4.0, 9.0, 5.0]
rowub = [ 5.0,  3.0,  9.0]
s = SCSSolver()
m = MathProgBase.ConicModel(s)
MathProgBase.loadproblem!(m, -obj, A, rowub, [(:NonNeg,1:3)],[(:NonNeg,1:5)])
MathProgBase.optimize!(m)
@test_approx_eq_eps MathProgBase.getobjval(m) -99.0 1e-3

# With eps = 1e-8, solution should be far more accurate
s = SCSSolver(eps=1e-8)
m = MathProgBase.ConicModel(s)
MathProgBase.loadproblem!(m, -obj, A, rowub, [(:NonNeg,1:3)],[(:NonNeg,1:5)])
MathProgBase.optimize!(m)
@test_approx_eq_eps MathProgBase.getobjval(m) -99.0 1e-5

# With a warmstart from the eps = 1e-8 solution, solution should be extremely accurate even after 1 iteration
push!(m.options, (:warm_start, true))
push!(m.options, (:max_iters, 1))
MathProgBase.optimize!(m)
@test_approx_eq_eps MathProgBase.getobjval(m) -99.0 1e-5

# Now let's do the same warmstart, but on a new instance of the same problem
primal_sol = m.primal_sol
dual_sol = m.dual_sol
slack = m.slack
s = SCSSolver(max_iters=1)
m = MathProgBase.ConicModel(s)
MathProgBase.loadproblem!(m, -obj, A, rowub, [(:NonNeg,1:3)],[(:NonNeg,1:5)])
MathProgBase.setwarmstart!(m, primal_sol; dual_sol = dual_sol, slack = slack)
MathProgBase.optimize!(m)
@test_approx_eq_eps MathProgBase.getobjval(m) -99.0 1e-5
