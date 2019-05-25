using Test
using SCS

tests = ["direct.jl",
         "indirect.jl",
         "options.jl"]

include("test_problems.jl")

for curtest in tests
    @testset "$curtest" begin
        include(curtest)
    end
end

include("MPB_wrapper.jl")
include("MOI_wrapper.jl")
