__precompile__()

module SCS

# To make Pkg aware that this dependency
# will be injected by BinaryProvider. 		
if VERSION > v"0.7.0-"
    using Libdl
end

if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    include("../deps/deps.jl")
else
    error("SCS not properly installed. Please run Pkg.build(\"SCS\") and restart julia")
end

using Compat

function __init__()
    vnum = VersionNumber(SCS_version())
    depsdir = realpath(joinpath(dirname(@__FILE__),"..","deps"))
    if vnum.major == 1
        error("Current SCS version installed is $(SCS_version()), but we require version 2.*.")
    end
end

include("types.jl")
include("c_wrapper.jl")
include("MPB_wrapper.jl")
include("MOI_wrapper.jl")

end # module
