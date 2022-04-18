if get(ENV, "BUILDKITE", "false") == "true"
    # This file requires a GPU in order to run. It gets tested as part of the
    # JuliaGPU CI on Buildkite. Contact @odow for more details.
    include("test_gpu.jl")
    exit(0)
end

using Test
using SCS

include("test_problems.jl")
@test SCS.scs_version() isa String
@test VersionNumber() >= v"3.2.0"
for solver in SCS.available_solvers
    @test SCS.scs_version(solver) isa String
    @test VersionNumber(SCS.scs_version(solver)) >= v"3.2.0"

    feasible_basic_problems(solver)
    test_options(solver)
end

include("MOI_wrapper.jl")
