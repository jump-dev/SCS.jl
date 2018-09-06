using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIB = MOI.Bridges

const MOIU = MOI.Utilities
MOIU.@model SCSModelData () (EqualTo, GreaterThan, LessThan) (Zeros, Nonnegatives, Nonpositives, SecondOrderCone, ExponentialCone, PositiveSemidefiniteConeTriangle) () (SingleVariable,) (ScalarAffineFunction,) (VectorOfVariables,) (VectorAffineFunction,)
const optimizer = MOIU.CachingOptimizer(SCSModelData{Float64}(), SCS.Optimizer())

# linear9test needs 1e-3 with SCS < 2.0 and 5e-1 with SCS 2.0
# linear2test needs 1e-4
const config = MOIT.TestConfig(atol=1e-4, rtol=1e-4)

@testset "Continuous linear problems" begin
    # AlmostSuccess for linear9 with SCS 2
    MOIT.contlineartest(MOIB.SplitInterval{Float64}(optimizer), config, ["linear9"])
end

@testset "Continuous conic problems" begin
    MOIT.contconictest(MOIB.RootDet{Float64}(MOIB.LogDet{Float64}(optimizer)), config, ["rsoc", "geomean", "psds", "rootdet", "logdets"])
end
