if haskey(ENV, "GITLAB_CI")
    # This file requires a GPU in order to run. It gets tested as part of the
    # JuliaGPU CI on Gitlab. Contact @odow for more details.
    include("test_gpu.jl")
    exit(0)
end

using Test
using SCS

solvers = SCS.available_solvers

include("test_problems.jl")

function feasible_basic_problems(solver)
    A = reshape([1.0],(1,1))
    solution = SCS_solve(solver, 1, 1, A, [1.0], [1.0], 1, 0, Int[], Int[], 0, 0, Float64[]);
    @test solution.ret_val == 1

    feasible_basic_conic(solver)
    feasible_exponential_conic(solver)
    feasible_sdp_conic(solver)
    feasible_pow_conic(solver)
end

for s in solvers
    feasible_basic_problems(s)
end

include("options.jl")

include("MOI_wrapper.jl")
include("MPB_wrapper.jl")
