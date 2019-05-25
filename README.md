# SCS

[![Build Status](https://travis-ci.org/JuliaOpt/SCS.jl.svg?branch=master)](https://travis-ci.org/JuliaOpt/SCS.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/yb4yfg4oryw7yten/branch/master?svg=true)](https://ci.appveyor.com/project/mlubin/scs-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/JuliaOpt/SCS.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaOpt/SCS.jl?branch=master)

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

To install custom built SCS binaries set the environmental variable `JULIA_SCS_LIBRARY_PATH` and call `Pkg.build("SCS")`. For instance, if the libraries are installed in `/opt/lib` just call
```julia
ENV["JULIA_SCS_LIBRARY_PATH"]="/opt/lib"
Pkg.build("SCS")
```
Note that your custom build binaries need to be compiled with the option `DLONG=1`. For instance, a minimal compilation script would be
```bash
$ cd <scs_dir>
$ make DLONG=1
$ julia
julia> ENV["JULIA_SCS_LIBRARY_PATH"]="<scs_dir>/out"
] build SCS
```
where `<scs_dir>` is SCS's source directory.

If you do not want BinaryProvider to download the default binaries on install set  `JULIA_SCS_LIBRARY_PATH`  before calling `Pkg.add("SCS")`.

To switch back to the default binaries clear `JULIA_SCS_LIBRARY_PATH` and call `Pkg.build("SCS")`.

## Usage

### MathProgBase wrapper
SCS implements the solver-independent [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) interface, and so can be used within modeling software like [Convex](https://github.com/JuliaOpt/Convex.jl) and [JuMP](https://github.com/JuliaOpt/JuMP.jl). The solver object is called `SCSSolver`.

### Options
All SCS solver options can be set through the direct interface(documented below) and through MathProgBase.
The list of options is defined the [`scs.h` header](https://github.com/cvxgrp/scs/blob/58e9af926fabc6674a9f488d4e9761a4f0fc451c/include/scs.h#L43).
To use these settings you can either pass them as keyword arguments to `SCS_solve` (high level interface) or as arguments to the `SCSSolver` constructor (MathProgBase interface), e.g.
```julia
# Direct
solution = SCS_solve(m, n, A, ..., psize; max_iters=10, verbose=0);
# MathProgBase (with Convex)
m = solve!(problem, SCSSolver(max_iters=10, verbose=0))
```

Moreover, You may select one of the linear solvers to be used by `SCSSolver` via `linear_solver` keyword. The options available are `SCS.Indirect` (the default) and `SCS.Direct`.

### High level wrapper

The file [`high_level_wrapper.jl`](https://github.com/JuliaOpt/SCS.jl/blob/master/src/high_level_wrapper.jl) is thoroughly commented. Here is the basic usage

We assume we are solving a problem of the form
```
minimize        c' * x
subject to      A * x + s = b
                s in K
```
where K is a product cone of

- zero cones,
- linear cones `{ x | x >= 0 }`,
- second-order cones `{ (t,x) | ||x||_2 <= t }`,
- semi-definite cones `{ X | X psd }`,
- exponential cones `{(x,y,z) | y e^(x/y) <= z, y>0 }`, and
- power cone `{(x,y,z) | x^a * y^(1-a) >= |z|, x>=0, y>=0}`.

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
- `p` is the array of power cone parameters
- `options` is a dictionary of options (see above).

The function is

```julia
function SCS_solve(m::Int, n::Int, A::SCSVecOrMatOrSparse, b::Array{Float64,},
    c::Array{Float64,}, f::Int, l::Int, q::Array{Int,}, qsize::Int, s::Array{Int,},
    ssize::Int, ep::Int, ed::Int, p::Array{Float64,}, psize::Int; options...)
```

and it returns an object of type Solution, which contains the following fields

```julia
type Solution
  x::Array{Float64, 1}
  y::Array{Float64, 1}
  s::Array{Float64, 1}
  status::ASCIIString
  ret_val::Int
  ...
```

Where `x` stores the optimal value of the primal variable, `y` stores the optimal value of the dual variable, `s` is the slack variable, `status` gives information such as `solved`, `primal infeasible`, etc.

### Low level wrapper

The low level wrapper directly calls SCS and is also thoroughly documented in [low_level_wrapper.jl](https://github.com/JuliaOpt/SCS.jl/blob/master/src/low_level_wrapper.jl). The low level wrapper performs the pointer manipulation necessary for the direct C call.

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
solve!(p, SCSSolver())
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

m = Model(solver=SCSSolver())
@variable(m, 0 <= take[items] <= 1)  # Define a variable for each item
@objective(m, Max, sum( values[item] * take[item] for item in items))
@constraint(m, sum( weight[item] * take[item] for item in items) <= 3)
solve(m)
println(getvalue(take))
# [Bronze] = 0.9999999496295456
# [  Gold] = 0.99999492720597
# [Silver] = 0.4666851698368782
```
