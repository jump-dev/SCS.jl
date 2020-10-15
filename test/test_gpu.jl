# This file requires a GPU in order to run. It gets tested as part of the
# JuliaGPU CI on Gitlab. Contact @odow for more details.

# CUDA_jll is not in our Project.toml, so we need to install it on GITLAB_CI.
import Pkg
Pkg.add("CUDA_jll")

using CUDA_jll  # CUDA_jll must be loaded _before_ SCS!
using SCS

using Test

using MathOptInterface
const MOI = MathOptInterface

const CONFIG = MOI.Test.TestConfig(atol=1e-5)

function _new_optimizer()
    return MOI.Bridges.full_bridge_optimizer(
        MOI.Utilities.CachingOptimizer(
            MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
            SCS.Optimizer(linear_solver = SCS.GpuIndirectSolver, eps = 1e-6),
        ),
        Float64,
    )
end

@test SCS.GpuIndirectSolver in SCS.available_solvers

@testset "contconictest" begin
    optimizer = _new_optimizer()
    MOI.Test.contlineartest(optimizer, CONFIG)
end

@testset "contconictest" begin
    optimizer = _new_optimizer()
    MOI.Test.contconictest(optimizer, CONFIG, ["rootdets", "logdets"])
end
