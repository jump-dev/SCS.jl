using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIB = MOI.Bridges

const MOIU = MOI.Utilities
MOIU.@model(SCSModelData,
            (),
            (MOI.EqualTo, MOI.GreaterThan, MOI.LessThan),
            (MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives, MOI.SecondOrderCone,
             MOI.ExponentialCone, MOI.PositiveSemidefiniteConeTriangle),
            (),
            (MOI.SingleVariable,),
            (MOI.ScalarAffineFunction,),
            (MOI.VectorOfVariables,),
            (MOI.VectorAffineFunction,))
for T in [SCS.Direct, SCS.Indirect]
    optimizer = MOIU.CachingOptimizer(SCSModelData{Float64}(),
                                      SCS.Optimizer(linear_solver=T, eps=1e-8,
                                                    verbose=0))

    @testset "SolverName" begin
        @test MOI.get(optimizer, MOI.SolverName()) == "SCS"
    end

    @testset "supports_allocate_load" begin
        @test MOIU.supports_allocate_load(optimizer.optimizer, false)
        @test !MOIU.supports_allocate_load(optimizer.optimizer, true)
    end

    config = MOIT.TestConfig(atol=1e-5)

    @testset "Unit" begin
        MOIT.unittest(MOIB.SplitInterval{Float64}(optimizer), config,
                      [# Quadratic functions are not supported
                       "solve_qcp_edge_cases", "solve_qp_edge_cases",
                       # Integer and ZeroOne sets are not supported
                       "solve_integer_edge_cases", "solve_objbound_edge_cases"])
    end

    @testset "Continuous linear problems with $T" begin
        MOIT.contlineartest(MOIB.SplitInterval{Float64}(optimizer), config)
    end

    @testset "Continuous conic problems with $T" begin
        MOIT.contconictest(MOIB.RootDet{Float64}(MOIB.LogDet{Float64}(MOIB.GeoMean{Float64}(MOIB.RSOC{Float64}(optimizer)))),
                           config, ["psds", "rootdets", "logdets"])
    end
end
