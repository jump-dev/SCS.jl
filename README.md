# SCS

[![Build Status](https://github.com/jump-dev/SCS.jl/workflows/CI/badge.svg?branch=master)](https://github.com/jump-dev/SCS.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/SCS.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/SCS.jl)

Julia wrapper for the [SCS](https://github.com/cvxgrp/scs) splitting cone
solver. SCS can solve linear programs, second-order cone programs, semidefinite
programs, exponential cone programs, and power cone programs.

## Installation

Install SCS.jl using the Julia package manager:
```julia
import Pkg; Pkg.add("SCS")
```

In addition to installing the SCS.jl pacakge, this will also download and
install the SCS binaries. (You do not need to install SCS separately.) If you
require a custom build of SCS, see the **Custom installation** instructions
below.

## Usage

### Use with Convex.jl

This example shows how we can model a simple knapsack problem with Convex and
use SCS to solve it.
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

### Use with JuMP

This example shows how we can model a simple knapsack problem with JuMP and use
SCS to solve it.
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

## Options

All SCS solver options can be set through `Convex.jl` or `MathOptInterface.jl`.

For example:
```julia
model = Model(optimizer_with_attributes(SCS.Optimizer, "max_iters" => 10))

# via MathOptInterface:
optimizer = SCS.Optimizer()
MOI.set(optimizer, MOI.RawParameter("max_iters"), 10)
MOI.set(optimizer, MOI.RawParameter("verbose"), 0)
```

Common options are:
 * `max_iters`: the maximum number of iterations to take
 * `verbose`: turn printing on (`1`) or off (`0`)
See the [`glbopts.h` header](https://github.com/cvxgrp/scs/blob/3aaa93c7aa04c7001df5e51b81f21b126dfa99b3/include/glbopts.h#L35)
for other options.

Select one of the linear solvers using the `linear_solver` option. The values
available are `SCS.DirectSolver` (the default) and `SCS.IndirectSolver`. A third
option for using a GPU is experimental, see the section below.

#### SCS on GPU

An experimental `SCS.GpuIndirectSolver` can be used by either providing the
appropriate libraries in a custom installation, or via the default binaries.
Not that `CUDA_jll` must be installed and loaded **before* `SCS`.
```julia
julia> import Pkg

julia> Pkg.add(Pkg.PackageSpec(name = "CUDA_jll", version = "10.1"))

julia> using CUDA_jll  # This must be called before `using SCS`.

julia> using SCS

julia> SCS.available_solvers
3-element Array{DataType,1}:
 SCS.DirectSolver
 SCS.IndirectSolver
 SCS.GpuIndirectSolver

julia> optimizer = SCS.Optimizer();

julia> MOI.set(
           optimizer,
           MOI.RawParameter("linear_solver"),
           SCS.GpuIndirectSolver,
       )
```

### Low-level wrapper

SCS provides a low-level interface to solve a problem directly, without
interfacing through MathOptInterface.

We assume we are solving a problem of the form
```
minimize        1/2 * x' * P * x + c' * x
subject to      A * x + s = b
                s in K
```
where `K` is a product cone of:
- zero cone,
- positive orthant `{ x | x ≥ 0 }`,
- box cone `{ (t,x) | t*l ≤ x ≤ t*u}`,
- second-order cone (SOC) `{ (t,x) | ||x||_2 ≤ t }`,
- semi-definite cone (SDC) `{ X | X is psd }`,
- exponential cone `{ (x,y,z) | y e^(x/y) ≤ z, y>0 }`,
- power cone `{ (x,y,z) | x^a * y^(1-a) ≥ |z|, x ≥ 0, y ≥ 0 }`, and
- dual power cone `{ (u,v,w) | (u/a)^a * (v/(1-a))^(1-a) ≥ |w|, u ≥ 0, v ≥ 0 }`.

To solve this problem with SCS, call `SCS.scs_solve`:
```julia
function scs_solve(
    linear_solver::Type{<:LinearSolver},
    m::Integer,
    n::Integer,
    A::AbstractMatrix,
    b::Vector{Float64},
    c::Vector{Float64},
    z::Integer,
    l::Integer,
    q::Vector{<:Integer},
    s::Vector{<:Integer},
    ep::Integer,
    ed::Integer,
    p::Vector{Float64},
    primal_sol::Vector{Float64}=zeros(n),
    dual_sol::Vector{Float64}=zeros(m),
    slack::Vector{Float64}=zeros(m);
    options...,
)
```
where:
- `A` is a matrix with `m` rows and `n` cols
- `P` is a positive semidefinite symmetric `n × n` matrix
- `b` is a vector of length `m`
- `c` is a vector of length `n`
- `z` is the number of primal zero / dual free cones, i.e., primal equality
  constraints
- `l` is the number of linear cones
- `bu` is the vector of upper bounds for the box cone
- `bl` is the vector of lower bounds for the box cone
- `q` is the array of SOCs sizes
- `s` is the array of SDCs sizes
- `ep` is the number of primal exponential cones
- `ed` is the number of dual exponential cones
- `p` is the array of power cone parameters (±1, with negative values for the
  dual cone)
- `primal_sol = zeros(n)` is a vector to warmstart the primal variables
- `dual_sol = zeros(m)` is a vector to warmstart the dual variables
- `slack = zeros(m)` is a vector to warmstart the slack variables
- `options` is a dictionary of options (see above).

Note: to successfully warmstart the solver all `primal_sol`, `dual_sol` and
`slack` must be provided **and** you must pass the `warm_start = true` option.

`scs_solve` returns an object of type `Solution`, which contains the following
fields:
```julia
mutable struct Solution{T}
    x::Vector{Float64}
    y::Vector{Float64}
    s::Vector{Float64}
    info::ScsInfo{T}
    ret_val::T
end
```
where `x` stores the optimal value of the primal variable, `y` stores the
optimal value of the dual variable, `s` is the slack variable, and `info`
contains various information about the solve step.

`SCS.raw_status(::ScsInfo)::String` describes the status, e.g. `"solved"`,
`"infeasible"`, `"unbounded"`, etc. For the precise return status, the value of
`ret_val` field should be compared with `https://github.com/cvxgrp/scs/blob/3aaa93c7aa04c7001df5e51b81f21b126dfa99b3/include/glbopts.h#L18`.

## Custom Installation

Custom build binaries are needed when using non-standard compile options, or
non-official julia binaries. Special caution is required during the compilation
of the `scs` libraries to ensure proper options and linking:

  * `libscsdir` and `libscsindir` need to be compiled with `DLONG=1`.
  * (optional) `libscsgpuindir` needs to be compiled with `DLONG=0`

All of these libraries should be linked against the OpenBLAS library which julia
uses. For the official julia binaries this can be achieved by e.g.
```bash
cd SCS_SOURCE_DIR
make purge
make USE_OPENMP=1 BLAS64=1 BLASSUFFIX=_64_ DLONG=1 BLASLDFLAGS="-L$JULIA_BLAS_PATH -lopenblas64_" out/libscsdir.so out/libscsindir.so
make clean
make USE_OPENMP=1 BLAS64=1 BLASSUFFIX=_64_ DLONG=0 BLASLDFLAGS="-L$JULIA_BLAS_PATH -lopenblas64_" out/libscsgpuindir.so
```
where
 * `SCS_SOURCE_DIR` is the main directory of the source of `scs`, and
 * `JULIA_BLAS_PATH` is the path to the directory containing BLAS library used
   by `julia`.
   - Before `julia-1.3`: `abspath(joinpath(Sys.BINDIR, "..", "lib", "julia"))`),
   - afterwards: the path to `BLAS` library artifact, e.g.
     `joinpath(OpenBLAS_jll.artifact_dir, "lib", "julia")`

To use custom built SCS binaries with `SCS.jl` set the environment variable
`JULIA_SCS_LIBRARY_PATH` to `SCS_SOURCE_DIR/opt` and build `SCS.jl`:
```julia
ENV["JULIA_SCS_LIBRARY_PATH"]="SCS_SOURCE_DIR/out"
using Pkg; Pkg.build("SCS")
```

To switch back to the default binaries remove `JULIA_SCS_LIBRARY_PATH` from
`ENV` and call `Pkg.build("SCS")` again.
