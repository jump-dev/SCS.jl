module SCS

if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    include("../deps/deps.jl")
else
    error("SCS not properly installed. Please run Pkg.build(\"SCS\") and restart julia")
end

if VERSION >= v"0.4.0-dev+3844"
    import Base.Libdl: RTLD_LAZY, RTLD_DEEPBIND, RTLD_GLOBAL, dlopen
end

function __init__()
    # default binaries need access to Julia's lapack symbols
    @unix_only dlopen(Base.liblapack_name, RTLD_LAZY|RTLD_DEEPBIND|RTLD_GLOBAL)
end

include("types.jl")
include("low_level_wrapper.jl")
include("high_level_wrapper.jl")
include("SCSSolverInterface.jl")  # MathProgBase interface

end # module
