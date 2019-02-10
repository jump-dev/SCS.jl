using Compat.Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

MOIU.@model(ModelData, (), (),
            (MOI.Zeros, MOI.Nonnegatives, MOI.SecondOrderCone,
             MOI.ExponentialCone, MOI.PositiveSemidefiniteConeTriangle),
            (), (), (), (MOI.VectorOfVariables,), (MOI.VectorAffineFunction,))

# UniversalFallback is needed for starting values
const cache = MOIU.UniversalFallback(ModelData{Float64}())

import SCS

for T in [SCS.Direct, SCS.Indirect]
    optimizer = SCS.Optimizer(linear_solver=T, eps=1e-8, verbose=0)
    MOI.empty!(cache)
    cached = MOIU.CachingOptimizer(cache, optimizer)

    # Essential bridges that are needed for all tests
    bridged = MOIB.Vectorize{Float64}(MOIB.NonposToNonneg{Float64}(cached))

    @testset "SolverName" begin
        @test MOI.get(optimizer, MOI.SolverName()) == "SCS"
    end

    @testset "supports_allocate_load" begin
        @test MOIU.supports_allocate_load(optimizer, false)
        @test !MOIU.supports_allocate_load(optimizer, true)
    end

    config = MOIT.TestConfig(atol=1e-5)

    @testset "Unit" begin
        MOIT.unittest(MOIB.SplitInterval{Float64}(bridged), config,
                      [# Quadratic functions are not supported
                       "solve_qcp_edge_cases", "solve_qp_edge_cases",
                       # Integer and ZeroOne sets are not supported
                       "solve_integer_edge_cases", "solve_objbound_edge_cases"])
    end

    @testset "Continuous linear problems with $T" begin
        MOIT.contlineartest(MOIB.SplitInterval{Float64}(bridged), config)
    end

    @testset "Continuous conic problems with $T" begin
        MOIT.contconictest(MOIB.RootDet{Float64}(MOIB.LogDet{Float64}(MOIB.GeoMean{Float64}(MOIB.RSOC{Float64}(bridged)))),
                           config, ["psds", "rootdets", "logdets"])
    end
end
