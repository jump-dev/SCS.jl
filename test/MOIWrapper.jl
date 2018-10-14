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

    @testset "Continuous linear problems with $T" begin
        MOIT.contlineartest(MOIB.SplitInterval{Float64}(optimizer), config)
    end

    @testset "Continuous conic problems with $T" begin
        MOIT.contconictest(MOIB.RootDet{Float64}(MOIB.LogDet{Float64}(optimizer)), config, ["rsoc", "geomean", "psds", "rootdet", "logdets"])
    end
end
