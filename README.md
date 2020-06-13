# SCS

[![Build Status](https://travis-ci.org/jump-dev/SCS.jl.svg?branch=master)](https://travis-ci.org/jump-dev/SCS.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/yb4yfg4oryw7yten/branch/master?svg=true)](https://ci.appveyor.com/project/mlubin/scs-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/jump-dev/SCS.jl/badge.svg?branch=master)](https://coveralls.io/r/jump-dev/SCS.jl?branch=master)

Julia wrapper for the [SCS](https://github.com/cvxgrp/scs) splitting cone
solver. SCS can solve linear programs, second-order cone programs, semidefinite
programs, exponential cone programs, and power cone programs.

## Installation

You can install SCS.jl through the Julia package manager:
```julia
julia> Pkg.add("SCS")
```

SCS.jl will use [BinaryProvider.jl](https://github.com/JuliaPackaging/BinaryProvider.jl) to automatically install the SCS binaries. Note that if you are not using the official Julia binaries from `https://julialang.org/downloads/` you may need a custom install of the SCS binaries.

## Custom Installation

Custom build binaries will allow to use e.g. the indirect solver on (a CUDA-enabled) gpu,
however special caution is required during the compilation of the `scs` libraries to ensure proper options and linking:

  * `libscsdir` and `libscsindir` need to be compiled with `DLONG=1`.
  * (optional) `libscsgpu` needs to be compiled with `DLONG=0`

All of these libraries should be linked against the OpenBLAS library which julia uses.
For the official julia binaries this can be achieved by e.g.

```bash
cd SCS_SOURCE_DIR
make purge
make USE_OPENMP=1 BLAS64=1 BLASSUFFIX=_64_ DLONG=1 BLASLDFLAGS="-L$JULIA_LIBRARY_PATH -lopenblas64_" out/libscsdir.so out/libscsindir.so
make clean
make USE_OPENMP=1 BLAS64=1 BLASSUFFIX=_64_ DLONG=0 BLASLDFLAGS="-L$JULIA_LIBRARY_PATH -lopenblas64_" out/libscsgpu.so
```
where
 * `SCS_SOURCE_DIR` is the main directory of the source of `scs`, and
 * `JULIA_LIBRARY_PATH` is the path to julia-shipped libraries (e.g. `abspath(joinpath(Sys.BINDIR, "..", "lib", "julia"))`)

To use custom built SCS binaries with `SCS.jl` set the environment variable
`JULIA_SCS_LIBRARY_PATH` to `SCS_SOURCE_DIR/opt` and build `SCS.jl`:
```julia
ENV["JULIA_SCS_LIBRARY_PATH"]="<scs_source_dir>/out"
using Pkg; Pkg.build("SCS")
```

To switch back to the default binaries delete `JULIA_SCS_LIBRARY_PATH` and call `Pkg.build("SCS")` again.

## Usage

### High-level interfaces
SCS implements the solver-independent [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) interface, and so can be used within modeling softwares like [Convex](https://github.com/JuliaOpt/Convex.jl) and [JuMP](https://github.com/jump-dev/JuMP.jl). The optimizer constructor is `SCS.Optimizer`.

A legacy [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) interface is available as well, in maintanence mode only.

### Options
All SCS solver options can be set through the direct interface(documented below), through `Convex.jl` or `MathOptInterface.jl`.
The list of options follows the [`glbopts.h` header](https://github.com/cvxgrp/scs/blob/0fd7ea85e8b0d878cacf5b1dbce557b330422ff7/include/glbopts.h#L30) in lowercase.
To use these settings you can either pass them as keyword arguments to `SCS_solve` (high level interface) or using the `SCS.Optimizer` constructor (MathOptInterface), e.g.
```julia
# Direct
solution = SCS_solve(m, n, A, ..., psize; max_iters=10, verbose=0);
# via MathOptInterface:
optimizer = SCS.Optimizer()
MOI.set(optimizer, MOI.RawParameter("max_iters"), 10)
MOI.set(optimizer, MOI.RawParameter("verbose"), 0)
```
or via specific helper functions:
```julia
problem = ... # JuMP problem
optimizer_constructor = optimizer_with_attributes(SCS.Optimizer, "max_iters" => 10, "verbose" => 0)
set_optimizer(problem, optimizer_constructor)
optimize!(problem)
```

Moreover, You may select one of the linear solvers to be used by `SCS.Optimizer` via `linear_solver` keyword.
The options available are `SCS.IndirectSolver` (the default) and `SCS.DirectSolver`.
An experimental `SCS.IndirectGpuSolver` can be used only with custom installation.

### High level wrapper

The file [`c_wrapper.jl`](https://github.com/jump-dev/SCS.jl/blob/master/src/c_wrapper.jl) is thoroughly commented. Here is the basic usage.

We assume we are solving a problem of the form
```
minimize        c' * x
subject to      A * x + s = b
                s in K
```
where `K` is a product cone of

- zero cones,
- positive orthant `{ x | x >= 0 }`,
- second-order cones (SOC) `{ (t,x) | ||x||_2 <= t }`,
- semi-definite cones (SDC) `{ X | X is psd }`,
- exponential cones `{ (x,y,z) | y e^(x/y) <= z, y>0 }`,
- power cone `{ (x,y,z) | x^a * y^(1-a) >= |z|, x>=0, y>=0 }`, and
- dual power cone `{ (u,v,w) | (u/a)^a * (v/(1-a))^(1-a) >= |w|, u>=0, v>=0 }`.

The problem data are

- `A` is the matrix with m rows and n cols
- `b` is of length m x 1
- `c` is of length n x 1
- `f` is the number of primal zero / dual free cones, i.e. primal equality constraints
- `l` is the number of linear cones
- `q` is the array of SOCs sizes
- `s` is the array of SDCs sizes
- `ep` is the number of primal exponential cones
- `ed` is the number of dual exponential cones
- `p` is the array of power cone parameters (±1, with negative values for the dual cone)
- `options` is a dictionary of options (see above).

The function is

```julia
function SCS_solve(linear_solver::Type{<:LinearSolver},
        m::Integer, n::Integer,
        A::SCS.VecOrMatOrSparse, b::Vector{Float64}, c::Vector{Float64},
        f::Integer, l::Integer, q::Vector{<:Integer}, s::Vector{<:Integer},
        ep::Integer, ed::Integer, p::Vector{Float64},
        primal_sol::Vector{Float64}=zeros(n),
        dual_sol::Vector{Float64}=zeros(m),
        slack::Vector{Float64}=zeros(m);
        options...)
```

and it returns an object of type `Solution`, which contains the following fields

```julia
mutable struct Solution{T<:SCSInt}
    x::Array{Float64, 1}
    y::Array{Float64, 1}
    s::Array{Float64, 1}
    info::SCSInfo{T}
    ret_val::T
end
```

Where `x` stores the optimal value of the primal variable, `y` stores the optimal value of the dual variable, `s` is the slack variable, and `info` contains various information about the solve step.
E.g. `SCS.raw_status(::SCSInfo)::String` describes the status, e.g. 'Solved', 'Intedeterminate', 'Infeasible/Inaccurate', etc.

### Convex and JuMP examples
This example shows how we can model a simple knapsack problem with Convex and use SCS to solve it.
```julia
using Convex, SCS
items  = [:Gold, :Silver, :Bronze]
values = [5.0, 3.0, 1.0]
weights = [2.0, 1.5, 0.3]

# Define a variable of size 3, each index representing an item
x = Variable(3)
p = maximize(x' * values, 0 <= x, x <= 1, x' * weights <= 3)
solve!(p, SCS.Optimizer)
println([items x.value])

# [:Gold 0.9999971880377178
#  :Silver 0.46667637765641057
#  :Bronze 0.9999998036351865]
```

This example shows how we can model a simple knapsack problem with JuMP and use SCS to solve it.
```julia
using JuMP, SCS
items  = [:Gold, :Silver, :Bronze]
values = Dict(:Gold => 5.0,  :Silver => 3.0,  :Bronze => 1.0)
weight = Dict(:Gold => 2.0,  :Silver => 1.5,  :Bronze => 0.3)

model = Model(SCS.Optimizer)
@variable(model, 0 <= take[items] <= 1)  # Define a variable for each item
@objective(model, Max, sum(values[item] * take[item] for item in items))
@constraint(model, sum(weight[item] * take[item] for item in items) <= 3)
optimize!(model)
println(value.(take))
# 1-dimensional DenseAxisArray{Float64,1,...} with index sets:
#     Dimension 1, Symbol[:Gold, :Silver, :Bronze]
# And data, a 3-element Array{Float64,1}:
#  1.0000002002226671
#  0.4666659513182934
#  1.0000007732744878

```
