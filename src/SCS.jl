module SCS

if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    include("../deps/deps.jl")
else
    error("SCS not properly installed. Please run Pkg.build(\"SCS\") and restart julia")
end

if VERSION >= v"0.4.0-dev+3844"
    import Base.Libdl: RTLD_LAZY, RTLD_DEEPBIND, RTLD_GLOBAL, dlopen
end

include("version.jl")

function __init__()
    # default binaries need access to Julia's lapack symbols
    @unix_only dlopen(Base.liblapack_name, RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
    if SCS_version() != scs_version
        error("Current SCS version installed is $(SCS_version()), but we require version $scs_version. On Linux and Windows, delete the contents of the `~/.julia/v0.3/SCS/deps` directory except for the files `build.jl` and `.gitignore`, then rerun Pkg.build(\"SCS\"). On OS X, run `using Homebrew; Homebrew.update()` in Julia.")
    end
end

include("types.jl")
include("low_level_wrapper.jl")
include("high_level_wrapper.jl")
include("SCSSolverInterface.jl")  # MathProgBase interface

end # module
