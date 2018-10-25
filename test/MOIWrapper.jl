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
    optimizer = MOIU.CachingOptimizer(SCSModelData{Float64}(), SCS.Optimizer(linear_solver=T, eps=1e-8))

    config = MOIT.TestConfig(atol=1e-5)

    @testset "Unit" begin
        MOIT.unittest(MOIB.SplitInterval{Float64}(optimizer),
                      # solve_blank_obj needs 1e-2 tolerance
                      MOIT.TestConfig(atol=1e-2, rtol=1e-2),
                      [# Quadratic functions are not supported
                       "solve_qcp_edge_cases", "solve_qp_edge_cases",
                       # Integer and ZeroOne sets are not supported
                       "solve_integer_edge_cases", "solve_objbound_edge_cases"])
    end

    @testset "Continuous linear problems" begin
        MOIT.unittest(MOIB.SplitInterval{Float64}(optimizer),
                      # solve_blank_obj needs 1e-2 tolerance
                      MOIT.TestConfig(atol=1e-2, rtol=1e-2),
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
