using Base.Test

tests = ["direct.jl",
         "indirect.jl",
         "mpb_linear.jl",
         "options.jl"]

include("test_problems.jl")

for curtest in tests
    @testset "$curtest" begin
        include(curtest)
    end
end

import SCS
include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
coniclineartest(SCS.SCSSolver(), duals=true, tol=1e-2)
conicSOCtest(SCS.SCSSolver(), duals=true, tol=1e-2)
conicEXPtest(SCS.SCSSolver(), duals=true, tol=1e-2)
conicSDPtest(SCS.SCSSolver(), duals=true, tol=1e-2)
