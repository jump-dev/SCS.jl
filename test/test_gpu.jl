# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# This file requires a GPU in order to run. It gets tested as part of the
# JuliaGPU CI on Gitlab. Contact @odow for more details.

# SCS_GPU_jll is not in our Project.toml, so we need to install it on GITLAB_CI.
import Pkg
Pkg.add("SCS_GPU_jll")

using SCS, SCS_GPU_jll
using Test

include("test_problems.jl")

@testset "GpuIndirectSolver" begin
    @test SCS.is_available(SCS.GpuIndirectSolver)
    version = SCS.scs_version(SCS.GpuIndirectSolver)
    @test VersionNumber(version) >= v"3.2.0"
    feasible_basic_problems(SCS.GpuIndirectSolver)
end
