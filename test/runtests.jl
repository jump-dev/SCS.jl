# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

if get(ENV, "BUILDKITE", "false") == "true"
    # This file requires a GPU in order to run. It gets tested as part of the
    # JuliaGPU CI on Buildkite. Contact @odow for more details.
    include("test_gpu.jl")
    exit(0)
end

@static if Sys.islinux() && Sys.ARCH == :x86_64
    import Pkg
    Pkg.add("SCS_MKL_jll")
    using SCS_MKL_jll
end

using Test
using SCS

solvers_to_test = Any[SCS.DirectSolver, SCS.IndirectSolver]
if SCS.is_available(SCS.MKLDirectSolver)
    push!(solvers_to_test, SCS.MKLDirectSolver)
end

include("test_problems.jl")

@testset "test-problems.jl" begin
    @test SCS.scs_version() isa String
    @test VersionNumber(SCS.scs_version()) >= v"3.2.0"
    for solver in solvers_to_test
        @test SCS.is_available(solver)
        @test SCS.scs_version(solver) isa String
        @test VersionNumber(SCS.scs_version(solver)) >= v"3.2.0"
        feasible_basic_problems(solver)
        test_options(solver)
    end
end

include("MOI_wrapper.jl")
