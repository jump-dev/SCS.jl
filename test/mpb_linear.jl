#############################################################################
# SCS.jl
# Wrapper around the SCS solver https://github.com/cvxgrp/scs
# See http://github.com/karanveerm/SCS.jl
#############################################################################
# test/mpb.jl
# Test the MathProgBase.jl interface for the SCS.jl solver wrapper
#############################################################################

using Base.Test
using MathProgBase
using SCS

solver = SCSSolver()
objtol = 1e-4
primaltol = 1e-4

# Stub some methods
MathProgBase.getreducedcosts(m::SCS.SCSMathProgModel) = nothing
MathProgBase.getconstrduals(m::SCS.SCSMathProgModel) = nothing

# min -x
# s.t. 2x + y <= 1.5
# x,y >= 0
# solution is (0.75,0) with objval -0.75
sol = linprog([-1,0],[2 1],'<',1.5,solver)
@test sol.status == :Optimal
@test_approx_eq_eps sol.objval -0.75 objtol
@test_approx_eq_eps norm(sol.sol - [0.75,0.0]) 0 primaltol

sol = linprog([-1,0],sparse([2 1]),'<',1.5,solver)
@test sol.status == :Optimal
@test_approx_eq_eps sol.objval -0.75 objtol
@test_approx_eq_eps norm(sol.sol - [0.75,0.0]) 0 primaltol

# test infeasible problem:
# min x
# s.t. 2x+y <= -1
# x,y >= 0
sol = linprog([1.0,0.0],[2.0 1.0],'<',-1.0,solver)
@test sol.status == :Infeasible

# test unbounded problem:
# min -x-y
# s.t. -x+2y <= 0
# x,y >= 0
sol = linprog([-1.,-1.],[-1. 2.],'<',[0.0],solver)
@test sol.status == :Unbounded

# More complicated through loadproblem!
A = [1.0 1.0 0.0 0.0 0.0;
     0.0 1.0 0.0 0.0 1.0;
     0.0 0.0 1.0 1.0 1.0]
collb = [0.0, 0.0, 0.0, 0.0, 0.0]
colub = [Inf, Inf, Inf, Inf, Inf]
obj   = [3.0, 4.0, 4.0, 9.0, 5.0]
sense = :Max
rowlb = [-Inf, -Inf, -Inf]
rowub = [ 5.0,  3.0,  9.0]
s = SCSSolver()
m = MathProgBase.model(s)
MathProgBase.loadproblem!(m, A, collb, colub, obj, rowlb, rowub, sense)
MathProgBase.optimize!(m)
@test_approx_eq_eps MathProgBase.getobjval(m) 99.0 1e-3
x = MathProgBase.getsolution(m)
@test_approx_eq_eps x[1] 2.0 1e-4
@test_approx_eq_eps x[2] 3.0 1e-4
@test_approx_eq_eps x[3] 0.0 1e-4
@test_approx_eq_eps x[4] 9.0 1e-4
@test_approx_eq_eps x[5] 0.0 1e-4
