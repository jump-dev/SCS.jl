module SCS

if isfile(joinpath(Pkg.dir("SCS"), "deps", "deps.jl"))
    include("../deps/deps.jl")
else
    error("SCS not properly installed. Please run Pkg.build(\"SCS\") and restart julia")
end

include("types.jl")
include("low_level_wrapper.jl")
include("high_level_wrapper.jl")
include("SCSSolverInterface.jl")  # MathProgBase interface

end # module
