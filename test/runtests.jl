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
conicSDPtest(SCS.SCSSolver(), duals=true, tol=1e-2)
