__precompile__()

module SCS

import Compat: unsafe_string

if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    include("../deps/deps.jl")
else
    error("SCS not properly installed. Please run Pkg.build(\"SCS\") and restart julia")
end

import Base.Libdl: RTLD_LAZY, RTLD_DEEPBIND, RTLD_GLOBAL, dlopen

function __init__()
    vnum = VersionNumber(SCS_version())
    # default binaries need access to Julia's lapack symbols
    @unix_only dlopen(Base.liblapack_name, RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
    if !(vnum.major == 1 && vnum.minor == 1)
        depsdir = realpath(joinpath(dirname(@__FILE__),"..","deps"))
        error("Current SCS version installed is $(SCS_version()), but we require version 1.1.*. On Linux and Windows, delete the contents of the `$depsdir` directory except for the files `build.jl` and `.gitignore`, then rerun Pkg.build(\"SCS\"). On OS X, run `using Homebrew; Homebrew.update()` in Julia.")
    end
end

include("types.jl")
include("low_level_wrapper.jl")
include("high_level_wrapper.jl")
include("SCSSolverInterface.jl")  # MathProgBase interface

end # module
