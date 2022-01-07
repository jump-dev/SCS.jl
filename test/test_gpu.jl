# This file requires a GPU in order to run. It gets tested as part of the
# JuliaGPU CI on Gitlab. Contact @odow for more details.

# CUDA_jll is not in our Project.toml, so we need to install it on GITLAB_CI.
import Pkg
Pkg.add(Pkg.PackageSpec(name = "CUDA_jll", version = "10.1"))

using CUDA_jll  # CUDA_jll must be loaded _before_ SCS!
using SCS
using Test

include("test_problems.jl")

@testset "Intro" begin
    @test SCS.GpuIndirectSolver in SCS.available_solvers
    @test(
          @ccall(SCS.gpuindirect.scs_sizeof_int()::Csize_t) ==
          sizeof(SCS.scsint_t(SCS.GpuIndirectSolver)),
    )
end

@testset "basic_feasible_problems" begin
    feasible_basic_problems(SCS.GpuIndirectSolver)
end

# TODO(odow): consider re-enabling these
# include("MOI_wrapper.jl")
# moi_tests(SCS.GpuIndirectSolver)
