@testset "Test the MathProgBase wrapper with linprog" begin
    include("mpb_linear.jl")
end

using Compat.Pkg: dir
@testset "MathProgBase" begin
    include(joinpath(dir("MathProgBase"),"test","conicinterface.jl"))
    coniclineartest(SCS.SCSSolver(), duals=true, tol=1e-2)
    conicSOCtest(SCS.SCSSolver(), duals=true, tol=1e-2)
    conicEXPtest(SCS.SCSSolver(), duals=true, tol=1e-2)
    conicSDPtest(SCS.SCSSolver(), duals=true, tol=1e-2)
end
