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
In addition to installing the `SCS.jl` package, this will also download and
install the SCS binaries. (You do not need to install SCS separately.)

To use a custom binary, read the [Custom solver binaries](https://jump.dev/JuMP.jl/stable/developers/custom_solver_binaries/)
section of the JuMP documentation.

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
MOI.set(optimizer, MOI.RawOptimizerAttribute("max_iters"), 10)
MOI.set(optimizer, MOI.RawOptimizerAttribute("verbose"), 0)
```

Common options are:
 * `max_iters`: the maximum number of iterations to take
 * `verbose`: turn printing on (`1`) or off (`0`)
See the [`glbopts.h` header](https://github.com/cvxgrp/scs/blob/3aaa93c7aa04c7001df5e51b81f21b126dfa99b3/include/glbopts.h#L35)
for other options.

### Linear solvers

`SCS` uses a linear solver internally, see
[this section](https://www.cvxgrp.org/scs/linear_solver/index.html#linear-system-solver)
of `SCS` documentation. `SCS.jl` ships with
* `SCS.DirectSolver` (sparse direct, the default) and
* `SCS.LinearSolver` (sparse indirect, by conjugate gradient)
enabled.

The find currently available linear solvers one can inspect `SCS.available_solvers`:
```julia
julia> using SCS

julia> SCS.available_solvers
2-element Vector{DataType}:
 SCS.DirectSolver
 SCS.IndirectSolver
```

To select the linear solver of choice
 * pass the `linear_solver` option to `JuMP.optimizer_with_attributes`, or to `MOI.OptimizerWithAttributes`;
 * specify the solver as the first argument when using `scs_solve` directly (see scetion Low-level wrapper below).

#### SCS with MKL Pardiso linear solver

To enable the MKL Pardiso (direct sparse) solver one needs to load `MKL_jll`
**before** `SCS`:

```julia
julia> import Pkg

julia> Pkg.add(Pkg.PackageSpec(name = "MKL_jll", version = "2022.2"))

julia> using MKL_jll    # This must be called before `using SCS`.

julia> using SCS

julia> SCS.available_solvers
3-element Vector{DataType}:
 SCS.DirectSolver
 SCS.IndirectSolver
 SCS.MKLDirectSolver
```

The `MKLDirectSolver` is available on `Linux x86_64` platform only.

#### SCS with Sparse GPU indirect solver (CUDA only)

> Note: as of version 1.0 the support for the GPU solver is broken (see
[this issue](https://github.com/jump-dev/SCS.jl/issues/245)).

To enable the indirect linear solver on gpu one needs to load `CUDA_jll`
**before** `SCS`:

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
```

The `GpuIndirectSolver` is available on `Linux x86_64` platform only.

### Low-level wrapper

SCS.jl provides a low-level interface to solve a problem directly, without
interfacing through MathOptInterface.

**This is an advanced interface with a risk of incorrect usage. For new users,
we recommend that you use the JuMP or Convex interfaces instead.**

SCS solves a problem of the form:
```
minimize        1/2 * x' * P * x + c' * x
subject to      A * x + s = b
                s in K
```
where `K` is a product cone of:
- zero cone
- positive orthant `{ x | x ≥ 0 }`
- box cone `{ (t,x) | t*l ≤ x ≤ t*u}`
- second-order cone (SOC) `{ (t,x) | ||x||_2 ≤ t }`
- semi-definite cone (SDC) `{ X | X is psd }`
- exponential cone `{ (x,y,z) | y e^(x/y) ≤ z, y>0 }`
- power cone `{ (x,y,z) | x^a * y^(1-a) ≥ |z|, x ≥ 0, y ≥ 0 }`
- dual power cone `{ (u,v,w) | (u/a)^a * (v/(1-a))^(1-a) ≥ |w|, u ≥ 0, v ≥ 0 }`.

To solve this problem with SCS, call `SCS.scs_solve`; see the docstring for
details.
