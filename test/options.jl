#############################################################################
# SCS.jl
# Wrapper around the SCS solver https://github.com/cvxgrp/scs
# See http://github.com/jump-dev/SCS.jl
#############################################################################
# test/options.jl
# Tests the ability to pass options
#############################################################################

using MathProgBase

@testset "SCS options" begin
    # The normal test
    A = [1.0 1.0 0.0 0.0 0.0;
        0.0 1.0 0.0 0.0 1.0;
        0.0 0.0 1.0 1.0 1.0]
    collb = [0.0, 0.0, 0.0, 0.0, 0.0]
    obj   = [3.0, 4.0, 4.0, 9.0, 5.0]
    rowub = [ 5.0,  3.0,  9.0]
    s = SCSSolver(eps_abs=1e-6, eps_rel=1e-7)
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m, -obj, A, rowub, [(:NonNeg,1:3)],[(:NonNeg,1:5)])
    MathProgBase.optimize!(m)

    @test isapprox(MathProgBase.getobjval(m), -99.0, atol=1e-5, rtol=0.0)
    @test !isapprox(MathProgBase.getobjval(m), -99.0, atol=1e-6, rtol=0.0)

    # With eps = 1e-10, solution should be far more accurate
    s = SCSSolver(eps_abs=1e-10, eps_rel=1e-12, acceleration_lookback=0)
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m, -obj, A, rowub, [(:NonNeg,1:3)],[(:NonNeg,1:5)])
    MathProgBase.optimize!(m)
    @test isapprox(MathProgBase.getobjval(m), -99.0, atol=1e-9, rtol=0.0)
    @test !isapprox(MathProgBase.getobjval(m), -99.0, atol=1e-10, rtol=0.0)

    # With a warmstart from the eps = 1e-10 solution, solution should be extremely accurate even after 1 iteration
    SCS.addoption!(m, :warm_start, true)
    SCS.addoption!(m, :max_iters, 151) #!!!! below 151 it gets Nan
    MathProgBase.optimize!(m)
    @test isapprox(MathProgBase.getobjval(m), -99.0, atol=1e-9, rtol=0.0)
    @test !isapprox(MathProgBase.getobjval(m), -99.0, atol=1e-10, rtol=0.0)

    # Now let's do the same warmstart, but on a new instance of the same problem
    primal_sol = m.primal_sol
    dual_sol = m.dual_sol
    slack = m.slack
    s = SCSSolver(max_iters=1)
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m, -obj, A, rowub, [(:NonNeg,1:3)],[(:NonNeg,1:5)])
    MathProgBase.setwarmstart!(m, primal_sol; dual_sol = dual_sol, slack = slack)
    MathProgBase.optimize!(m)
    @test isapprox(MathProgBase.getobjval(m), -99.0, atol=1e-9, rtol=0.0)
    @test !isapprox(MathProgBase.getobjval(m), -99.0, atol=1e-10, rtol=0.0)

    # tests for incorrect options
    s = SCSSolver(eps_abs=1e-12, epps=1.0)
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m, -obj, A, rowub, [(:NonNeg,1:3)],[(:NonNeg,1:5)])

    @test_throws ArgumentError MathProgBase.optimize!(m)

    err = try
        MathProgBase.optimize!(m)
    catch ex
        ex
    end
    @test err.msg == """Unrecognized option passed to SCS: epps;\nRecognized options are: linear_solver, normalize, scale, rho_x, max_iters, eps_abs, eps_rel, eps_infeas, alpha, time_limit_secs, verbose, warm_start, acceleration_lookback, acceleration_interval, adaptive_scaling, write_data_filename and log_csv_filename."""

    # tests for incorrect options
    s = SCSSolver(linear_solver="AAA", eps_abs=1e-12, epps=1.0)
    m = MathProgBase.ConicModel(s)
    MathProgBase.loadproblem!(m, -obj, A, rowub, [(:NonNeg,1:3)],[(:NonNeg,1:5)])

    @test_throws ArgumentError MathProgBase.optimize!(m)

    err = try
        MathProgBase.optimize!(m)
    catch ex
        ex
    end

    let msg = "Unrecognized linear_solver passed to SCS: AAA;\nRecognized options are: "
        if isdefined(SCS, :gpuindirect)
            msg *= "SCS.DirectSolver, SCS.IndirectSolver and SCS.GpuIndirectSolver."
        else
            msg *= "SCS.DirectSolver and SCS.IndirectSolver."
        end
        @test err.msg == msg
    end
end
