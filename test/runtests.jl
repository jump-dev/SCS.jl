tests = ["direct.jl",
         "mpb_linear.jl",
         "options.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end

import SCS
include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
coniclineartest(SCS.SCSSolver(), duals=true, tol=1e-2)
conicSOCtest(SCS.SCSSolver(), duals=true, tol=1e-2)
conicEXPtest(SCS.SCSSolver(), duals=true, tol=1e-2)
# TODO: duals don't work for SDPs
# TODO: this test fails, since the semidefinite cone sizes are wrong
#conicSDPtest(SCS.SCSSolver(), duals=false, tol=1e-2)
