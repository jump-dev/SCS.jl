module SCS

using Libdl

if haskey(ENV, "JULIA_SCS_LIBRARY_PATH") || VERSION < v"1.3"
    if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
        include(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    else
        error("SCS not properly installed. Please run Pkg.build(\"SCS\") and restart julia")
    end

    function __init__()
        vnum = VersionNumber(SCS_version())
        depsdir = realpath(joinpath(dirname(@__FILE__),"..","deps"))
        if vnum.major == 1 || (vnum.major == 2 && vnum.minor != 1)
            error("Current SCS version installed is $(SCS_version()), but we require version 2.1.*")
        end
    end
else
    import SCS_jll
    const indirect = SCS_jll.libscsindir
    const direct = SCS_jll.libscsdir
    if Sys.islinux() && Libdl.dlopen("libcuda.$(Libdl.dlext)"; throw_error=false) !== nothing
        import SCS_GPU_jll
	const gpuindirect = SCS_jll.libscsgpuindir
    end
end

include("types.jl")

const available_solvers = [DirectSolver, IndirectSolver]

include("c_wrapper.jl")
include("MPB_wrapper.jl")
include("MOI_wrapper.jl")

end # module
