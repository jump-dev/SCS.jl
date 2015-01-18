#############################################################################
# SCS.jl
# Wrapper around the SCS solver https://github.com/cvxgrp/scs
# See http://github.com/JuliaOpt/SCS.jl
#############################################################################
# test/mpb.jl
# Test the MathProgBase.jl interface for the SCS.jl solver wrapper
#############################################################################

using Base.Test
using MathProgBase.SolverInterface
using SCS

s = SCSSolver()
# Problem 1 - all vars in nonneg cone
# min -3x - 2y - 4z
# st    x +  y +  z == 3
#            y +  z == 2
#       x>=0 y>=0 z>=0
# Opt solution = -11
# x = 1, y = 0, z = 2
m = MathProgBase.model(s)
MathProgBase.loadconicproblem!(m,
                    [-3.0, -2.0, -4.0],
                    [ 1.0   1.0   1.0;
                      0.0   1.0   1.0],
                    [ 3.0,  2.0],
                    [(:Zero, 1:2)],
                    [(:NonNeg, 1:3)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) -11 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 0.0 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 2.0 1e-4

# Problem 2 - mixed free, nonneg, nonpos, zero, shuffled cones
# min  3x + 2y - 4z + 0s
# st    x           -  s  == -4    (i.e. x >= -4)
#            y            == -3
#       x      +  z       == 12
#       x free
#       y <= 0
#       z >= 0
#       s zero
# Opt solution = -82
# x = -4, y = -3, z = 16, s == 0
m = MathProgBase.model(s)
MathProgBase.loadconicproblem!(m,
                    [ 3.0,  2.0, -4.0,  0.0],
                    [ 1.0   0.0   0.0  -1.0;
                      0.0   1.0   0.0   0.0;
                      1.0   0.0   1.0   0.0],
                    [-4.0, -3.0, 12.0],
                    [(:Zero,1:3)],
                    [(:Free,1), (:NonNeg,3), (:Zero,4), (:NonPos,2)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) -82 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[1] -4.0 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[2] -3.0 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 16.0 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[4]  0.0 1e-4

# Problem 3 - SOC
# min 0x - 1y - 1z
#  st  x            == 1
#      x >= ||(y,z)||
m = MathProgBase.model(s)
MathProgBase.loadconicproblem!(m,
                    [ 0.0, -1.0, -1.0],
                    [ 1.0   0.0   0.0],
                    [ 1.0],
                    [(:Zero,1)],
                    [(:SOC,1:3)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) -sqrt(2.0) 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 1.0/sqrt(2.0) 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 1.0/sqrt(2.0) 1e-4

# Problem 4 - ExpPrimal
# min x + y + z
#  st  y e^(x/y) <= x, y > 0 (i.e (x, y, z) are in the exponential primal cone)
#      x == 1
#      y == 2
m = MathProgBase.model(s)
MathProgBase.loadconicproblem!(m,
[1.0, 1.0, 1.0],
[0.0 1.0 0.0;
 1.0 0.0 0.0],
[2.0, 1.0],
[(:Zero,1:2)],
[(:ExpPrimal, 1:3)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) 6.2974 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 2.0 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 3.29744 1e-4


function is_symmetric(A::Matrix)
  return all(A - A' .< 1e-4)
end

# Problem 5 - SDP
# min y[1, 2]
#  st y[2, 1] == 1
#     y in SDP cone
# If symmetricity constraint is working, y[1, 2] will be 1 else unbounded
m = MathProgBase.model(s)
c = [0, 1, 0, 0, 0, 0, 0, 0, 0];
A = -eye(9)
A = [A; [0, 0, 0, 1, 0, 0, 0, 0, 0]']
b = zeros(size(A, 1), 1)
b[10] = 1
MathProgBase.loadconicproblem!(m, c, A, b, [(:SDP, 1:9), (:Zero, 10)], [(:Free, 1:9)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test is_symmetric(reshape(MathProgBase.getsolution(m), 3, 3))

# Problem 6 - SDP
# Same as problem 5, except we enforce :SDP on the var_cone
m = MathProgBase.model(s)
c = [0, 1, 0, 0, 0, 0, 0, 0, 0];
A = [0, 0, 0, 1, 0, 0, 0, 0, 0]'
b = [1]
MathProgBase.loadconicproblem!(m, c, A, b, [(:Zero, 1:1)], [(:SDP, 1:9)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test is_symmetric(reshape(MathProgBase.getsolution(m), 3, 3))
