#############################################################################
# SCS.jl
# Wrapper around the SCS solver https://github.com/cvxgrp/scs
# See http://github.com/JuliaOpt/SCS.jl
#############################################################################
# test/mpb.jl
# Test the MathProgBase.jl interface for the SCS.jl solver wrapper
#############################################################################

using Base.Test
importall MathProgBase
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
                    [(:SOC,1:3)])
MathProgBase.optimize!(m)
@test MathProgBase.status(m) == :Optimal
@test_approx_eq_eps MathProgBase.getobjval(m) -sqrt(2.0) 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[1] 1.0 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[2] 1.0/sqrt(2.0) 1e-4
@test_approx_eq_eps MathProgBase.getsolution(m)[3] 1.0/sqrt(2.0) 1e-4
