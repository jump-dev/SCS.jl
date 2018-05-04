using Base.Test

tests = ["direct.jl",
         "indirect.jl",
         "options.jl"]

include("test_problems.jl")

for curtest in tests
    @testset "$curtest" begin
        include(curtest)
    end
end

include("MPBWrapper.jl")
include("MOIWrapper.jl")
