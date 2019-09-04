module SCS

if isfile(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
    include(joinpath(dirname(@__FILE__), "..", "deps", "deps.jl"))
else
    error("SCS not properly installed. Please run Pkg.build(\"SCS\") and restart julia")
end

function __init__()
    vnum = VersionNumber(SCS_version())
    depsdir = realpath(joinpath(dirname(@__FILE__), "..", "deps"))
    if vnum.major != 2 || vnum.minor != 0
        error("Current SCS version installed is $(SCS_version()), but we require version 2.0.*.")
    end
end

include("types.jl")
include("c_wrapper.jl")
include("MPB_wrapper.jl")
include("MOI_wrapper.jl")

end # module
