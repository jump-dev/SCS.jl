using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

# UniversalFallback is needed for starting values
const cache = MOIU.UniversalFallback(MOIU.Model{Float64}())

import SCS

for T in [SCS.Direct, SCS.Indirect]
    optimizer = SCS.Optimizer(linear_solver=T, eps=1e-8)
    MOI.set(optimizer, MOI.Silent(), true)

    @testset "SolverName" begin
        @test MOI.get(optimizer, MOI.SolverName()) == "SCS"
    end

    @testset "supports_allocate_load" begin
        @test MOIU.supports_allocate_load(optimizer, false)
        @test !MOIU.supports_allocate_load(optimizer, true)
    end

    MOI.empty!(cache)
    cached = MOIU.CachingOptimizer(cache, optimizer)
    bridged = MOIB.full_bridge_optimizer(cached, Float64)
    config = MOIT.TestConfig(atol=1e-5)

    @testset "Unit" begin
        MOIT.unittest(bridged, config, [
            # `TimeLimitSec` not supported.
            "time_limit_sec",
            # ArgumentError: The number of constraints in SCSModel must be greater than 0
            "solve_unbounded_model",
            # Quadratic functions are not supported
            "solve_qp_edge_cases",
            # Integer and ZeroOne sets are not supported
            "solve_integer_edge_cases", "solve_objbound_edge_cases",
            "solve_zero_one_with_bounds_1",
            "solve_zero_one_with_bounds_2",
            "solve_zero_one_with_bounds_3"])
    end

    @testset "Continuous linear problems with $T" begin
        MOIT.contlineartest(bridged, config)
    end

    @testset "Continuous quadratic problems with $T" begin
        MOIT.qcptest(bridged, config)
    end

    @testset "Continuous conic problems with $T" begin
        MOIT.contconictest(bridged, config, ["rootdets", "logdets"])
    end
end
