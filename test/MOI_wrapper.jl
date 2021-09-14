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
            # Unexpected test failures:
            #   Numerical rounding issues in returned function
            "test_basic_ScalarQuadraticFunction_LessThan",
            "test_basic_VectorOfVariables_GeometricMeanCone",
            #   UnsupportedAttribute not thrown
            "test_model_copy_to_UnsupportedAttribute",
            #   ConstraintDualStart not supported correctly
            "test_model_ModelFilter_AbstractConstraintAttribute",
            #   cone dimensions 241 not equal to num rows in A = m = 31
            "test_conic_LogDetConeTriangle",
            "test_conic_NormNuclearCone",
            "test_conic_NormSpectralCone",
            "test_conic_PositiveSemidefiniteConeSquare_3",
            "test_conic_PositiveSemidefiniteConeSquare_VectorAffineFunction",
            "test_conic_PositiveSemidefiniteConeSquare_VectorOfVariables",
            "test_conic_PositiveSemidefiniteConeTriangle",
            "test_conic_RootDetConeTriangle",
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
            # Upstream failures
            #   #1431
            "test_model_LowerBoundAlreadySet",
            "test_model_UpperBoundAlreadySet",
            #   #1602
            "test_model_ListOfConstraintAttributesSet",
            "test_objective_FEASIBILITY_SENSE_clears_objective",
            "test_conic_SecondOrderCone_negative_initial_bound",
            #   #1603
            "test_modification_set_function_single_variable",
            "test_objective_incorrect_modifications",
            #   #1604
            "test_solve_optimize_twice",
            #   #1605
            "test_basic_VectorAffineFunction_LogDetConeTriangle",
            "test_basic_VectorAffineFunction_RootDetConeTriangle",
            "test_basic_VectorOfVariables_LogDetConeTriangle",
            "test_basic_VectorOfVariables_RootDetConeTriangle",
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

end  # module

TestSCS.runtests()
