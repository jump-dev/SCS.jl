# SCS

[![Build Status](https://travis-ci.org/karanveerm/SCS.jl.png)](https://travis-ci.org/karanveerm/SCS.jl)
[![Coverage Status](https://coveralls.io/repos/karanveerm/SCS.jl/badge.png?branch=dev)](https://coveralls.io/r/karanveerm/SCS.jl?branch=dev)

Julia wrapper around [SCS](https://github.com/cvxgrp/scs) for [CVX.jl](https://github.com/cvxgrp/CVX.jl).

## High Level Wrapper

The file [`high_level_wrapper.jl`](https://github.com/karanveerm/SCS.jl/blob/master/src/high_level_wrapper.jl) is thoroughly commented. Here is the basic usage

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
- semi-definite cones `{ X | X psd }`, and
- exponential cones `{(x,y,z) | y e^(x/y) <= z, y>0 }`.

The other problem data are

- `A` is the matrix with m rows and n cols
- `b` is of length m x 1
- `c` is of length n x 1
- `f` is the number of primal zero / dual free cones, i.e. primal equality constraints
- `l` is the number of linear cones
- `q` is the array of SOCs sizes
- `s` is the array of SDCs sizes
- `ep` is the number of primal exponential cones
- `ed` is the number of dual exponential cones.

The function is

```
function SCS_solve(m::Int64, n::Int64, A::SCSVecOrMatOrSparse, b::Array{Float64,},
    c::Array{Float64,}, f::Clong, l::Clong, q::Array{Int64,}, qsize::Clong, s::Array{Int64,},
    ssize::Clong, ep::Clong, ed::Clong)
```

and it returns an object of type Solution, which contains the following fields

```
type Solution
  x::Array{Float64, 1}
  y::Array{Float64, 1}
  s::Array{Float64, 1}
  status::ASCIIString
  ret_val::Int64
  ...
```

Where `x` stores the optimal value of the primal variable, `y` stores the optimal value of the dual variable, `s` is the slack variable, `status` gives information such as `solved`, `primal infeasible`, etc.

## Low Level Wrapper

The low level wrapper directly calls SCS and is also thoroughly documented in [low_level_wrapper.jl](https://github.com/karanveerm/SCS.jl/blob/master/src/low_level_wrapper.jl). The low level wrapper performs the pointer manipulation necessary for the direct C call.
