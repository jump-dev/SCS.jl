module SCS

import MathOptInterface
import Requires
import SCS_jll
import SparseArrays

function __init__()
    global indirect = SCS_jll.libscsindir
    global direct = SCS_jll.libscsdir

    Requires.@require(
        CUDA_jll = "e9e359dc-d701-5aa8-82ae-09bbf812ea83",
        begin
            import SCS_GPU_jll
            global gpuindirect = SCS_GPU_jll.libscsgpuindir
            push!(available_solvers, GpuIndirectSolver)
        end
    )
    Requires.@require(
        MKL_jll = "856f044c-d86e-5d09-b602-aeab76dc8ba7",
        begin
            import SCS_MKL_jll
            import SCS_MKL_jll.MKL_jll
            import Libdl
            global mkldirect = SCS_MKL_jll.libscsmkl

            MKLROOT = joinpath(MKL_jll.artifact_dir, "lib")
            # dlopen with RTLD_GLOBAL to avoid undefined symbols when dlopening
            # `lib_mkl_gnu_thread.so`
            Libdl.dlopen(
                joinpath(MKLROOT, "libmkl_core.so"),
                Libdl.RTLD_LAZY | Libdl.RTLD_DEEPBIND | Libdl.RTLD_GLOBAL,
            )
            Libdl.dlopen(joinpath(MKLROOT, "libmkl_gnu_thread.so"))

            push!(available_solvers, MKLDirectSolver)
        end
    )
    return
end

include("c_wrapper.jl")
include("linear_solvers/direct.jl")
include("linear_solvers/indirect.jl")
include("linear_solvers/gpu_indirect.jl")
include("linear_solvers/mkl_direct.jl")
include("MOI_wrapper/MOI_wrapper.jl")

const available_solvers = [DirectSolver, IndirectSolver]

export scs_solve

end
