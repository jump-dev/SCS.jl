tests = ["direct.jl",
         "mpb_linear.jl",
         "mpb_conic.jl",
         "options.jl"]

println("Running tests:")

for curtest in tests
    println(" Test: $(curtest)")
    include(curtest)
end
