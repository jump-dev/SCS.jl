__precompile__()

module SCS

if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    include("../deps/deps.jl")
else
    error("SCS not properly installed. Please run Pkg.build(\"SCS\") and restart julia")
end

import Base.Libdl: RTLD_LAZY, RTLD_DEEPBIND, RTLD_GLOBAL, dlopen

function __init__()
    vnum = VersionNumber(SCS_version())
    # default binaries need access to Julia's lapack symbols
    if is_unix()
        dlopen(Base.liblapack_name, RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
    end
    depsdir = realpath(joinpath(dirname(@__FILE__),"..","deps"))
    if vnum.major == 1
        error("Current SCS version installed is $(SCS_version()), but we require version 2.*. On Linux and Windows, delete the contents of the `$depsdir` directory except for the files `build.jl` and `.gitignore`, then rerun Pkg.build(\"SCS\"). On OS X, run `using Homebrew; Homebrew.update()` in Julia.")
    end
end

include("types.jl")
include("c_wrapper.jl")
include("SCSSolverInterface.jl")  # MathProgBase interface

end # module
