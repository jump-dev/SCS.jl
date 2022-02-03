if get(ENV, "BUILDKITE", "false") == "true"
    # This file requires a GPU in order to run. It gets tested as part of the
    # JuliaGPU CI on Buildkite. Contact @odow for more details.
    include("test_gpu.jl")
    exit(0)
end

using Test
using SCS

include("test_problems.jl")
for s in SCS.available_solvers
    @test @ccall(SCS.direct.scs_sizeof_int()::Csize_t) ==
          sizeof(SCS.scsint_t(s))
    feasible_basic_problems(s)
    test_options(s)
end

include("MOI_wrapper.jl")
