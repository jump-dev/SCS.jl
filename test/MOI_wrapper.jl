module TestSCS

using Test
using MathOptInterface
import SCS

const MOI = MathOptInterface

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

test_DirectSolver() = _test_runtests(SCS.DirectSolver)

test_IndirectSolver() = _test_runtests(SCS.IndirectSolver)

function _test_runtests(linear_solver)
    optimizer = SCS.Optimizer()
    MOI.set(
        optimizer,
        MOI.RawOptimizerAttribute("linear_solver"),
        linear_solver,
    )
    MOI.set(optimizer, MOI.RawOptimizerAttribute("eps"), 1e-6)
    MOI.set(optimizer, MOI.Silent(), true)
    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            optimizer,
        ),
        Float64,
    )
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            atol = 1e-4,
            exclude = Any[
                MOI.ConstraintBasisStatus,
                MOI.VariableBasisStatus,
                MOI.ConstraintName,
                MOI.VariableName,
                MOI.ObjectiveBound,
            ],
        ),
        exclude = String[
            # TODO(odow): unimplemented features
            "test_attribute_SolverVersion",
            # Expected test failures:
            #   ArgumentError: The number of constraints must be greater than 0
            "test_attribute_RawStatusString",
            "test_attribute_SolveTimeSec",
            "test_objective_ObjectiveFunction_blank",
            "test_solve_TerminationStatus_DUAL_INFEASIBLE",
            #   Problem is a nonconvex QP
            "test_basic_ScalarQuadraticFunction_EqualTo",
            "test_basic_ScalarQuadraticFunction_GreaterThan",
            "test_basic_ScalarQuadraticFunction_Interval",
            "test_basic_VectorQuadraticFunction_",
            "test_quadratic_SecondOrderCone_basic",
            "test_quadratic_nonconvex_",
            #   power cone error, values must be in [-1,1]
            "test_conic_DualPowerCone_VectorOfVariables",
            "test_conic_DualPowerCone_VectorAffineFunction",
            "test_conic_PowerCone_VectorAffineFunction",
            "test_conic_PowerCone_VectorOfVariables",
            #   MathOptInterface.jl issue #1431
            "test_model_LowerBoundAlreadySet",
            "test_model_UpperBoundAlreadySet",
        ],
    )
    return
end

function test_RawOptimizerAttribute()
    model = SCS.Optimizer()
    MOI.set(model, MOI.RawOptimizerAttribute("eps"), 1.0)
    @test MOI.get(model, MOI.RawOptimizerAttribute("eps")) == 1.0
    @test MOI.get(model, MOI.RawOptimizerAttribute("eps")) == 1.0
    MOI.set(model, MOI.RawOptimizerAttribute("eps"), 2.0)
    @test MOI.get(model, MOI.RawOptimizerAttribute("eps")) == 2.0
    @test MOI.get(model, MOI.RawOptimizerAttribute("eps")) == 2.0
end

function test_constrained_variables()
    model = MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            SCS.Optimizer(),
        ),
        Float64,
    )
    @test MOI.supports_constraint(
        model,
        MOI.VectorOfVariables,
        MOI.PositiveSemidefiniteConeTriangle,
    )
    x = MOI.add_variables(model, 6)
    f = MOI.VectorOfVariables(x)
    s = MOI.PositiveSemidefiniteConeTriangle(3)
    @test isa(
        MOI.add_constraint(model, f, s),
        MOI.ConstraintIndex{typeof(f),typeof(s)},
    )
    return
end

function test_unsupported()
    model = MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}())
    optimizer = SCS.Optimizer()
    x = MOI.add_variable(model)
    MOI.add_constraint(model, 1.0x, MOI.EqualTo(1.0))
    err = MOI.UnsupportedConstraint{
        MOI.ScalarAffineFunction{Float64},
        MOI.EqualTo{Float64},
    }()
    @test_throws err MOI.optimize!(optimizer, model)
    MOI.empty!(model)
    x = MOI.add_variable(model)
    MOI.set(model, MOI.Test.UnknownVariableAttribute(), x, 1.0)
    err = MOI.UnsupportedAttribute{MOI.Test.UnknownVariableAttribute}
    @test_throws err MOI.optimize!(optimizer, model)
    MOI.empty!(model)
    x = MOI.add_variable(model)
    c = MOI.add_constraint(model, MOI.Utilities.vectorize([1.0x]), MOI.Zeros(1))
    MOI.set(model, MOI.Test.UnknownConstraintAttribute(), c, 1.0)
    err = MOI.UnsupportedAttribute{MOI.Test.UnknownConstraintAttribute}
    @test_throws err MOI.optimize!(optimizer, model)
end

end  # module

TestSCS.runtests()
