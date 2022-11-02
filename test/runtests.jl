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

using Test
# MKL_jll is not in our Project.toml, so we need to install it.
import Pkg
Pkg.add(Pkg.PackageSpec(name = "MKL_jll", version = "2022"))
using MKL_jll  # MKL_jll must be loaded _before_ SCS!
using SCS

if Sys.islinux() && Sys.ARCH == :x86_64
    @test SCS.MKLDirectSolver in SCS.available_solvers
end

include("test_problems.jl")
@test SCS.scs_version() isa String
@test VersionNumber(SCS.scs_version()) >= v"3.2.0"
for solver in SCS.available_solvers
    @test SCS.scs_version(solver) isa String
    @test VersionNumber(SCS.scs_version(solver)) >= v"3.2.0"

    feasible_basic_problems(solver)
    test_options(solver)
end

include("MOI_wrapper.jl")
