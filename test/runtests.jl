import Pkg
Pkg.add(Pkg.PackageSpec(name = "MathOptInterface", rev = "master"))

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
    feasible_basic_problems(s)
end

include("MOI_wrapper.jl")
