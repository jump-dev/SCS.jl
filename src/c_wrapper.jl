# Copyright (c) 2014: SCS.jl contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

abstract type AbstractSCSType end

mutable struct ScsMatrix{T} <: AbstractSCSType
    values::Ptr{Cdouble}
    rowval::Ptr{T}
    colptr::Ptr{T}
    m::T
    n::T
end

mutable struct ScsSettings{T} <: AbstractSCSType
    normalize::T # boolean, heuristic data rescaling
    scale::Cdouble # if normalized, rescales by this factor
    adaptive_scale::T # boolean, whether to adaptively update `scale`
    rho_x::Cdouble # x equality constraint scaling
    max_iters::T # maximum iterations to take
    eps_abs::Cdouble # absolute convergence tolerance
    eps_rel::Cdouble # relative convergence tolerance
    eps_infeas::Cdouble # infeasible convergence tolerance
    alpha::Cdouble # Douglas-Rachford relaxation parameter
    time_limit_secs::Cdouble # time limit in seconds
    verbose::T # boolean, write out progress
    warm_start::T # boolean, warm start (put initial guess in Sol struct)
    acceleration_lookback::T # memory for acceleration
    acceleration_interval::T # interval to apply acceleration
    write_data_filename::Cstring # if set dump raw problem data to the file
    log_csv_filename::Cstring # if set log solve data to the file
    function ScsSettings(linear_solver::Type{<:LinearSolver})
        settings = new{scsint_t(linear_solver)}()
        scs_set_default_settings(linear_solver, settings)
        return settings
    end
end

mutable struct ScsData{T} <: AbstractSCSType
    m::T
    n::T
    A::Ptr{ScsMatrix{T}}
    P::Ptr{Cvoid} # Ptr{ScsMatrix{T}} or C_NULL
    b::Ptr{Cdouble}
    c::Ptr{Cdouble}
end

mutable struct ScsCone{T} <: AbstractSCSType
    z::T
    l::T
    bu::Ptr{Cdouble}
    bl::Ptr{Cdouble}
    bsize::T
    q::Ptr{T}
    qsize::T
    s::Ptr{T}
    ssize::T
    cs::Ptr{T}
    cssize::T
    ep::T
    ed::T
    p::Ptr{Cdouble}
    psize::T
    d::Ptr{T}
    dsize::T
    nuc_m::Ptr{T}
    nuc_n::Ptr{T}
    nucsize::T
    ell1::Ptr{T}
    ell1_size::T
    sl_n::Ptr{T}
    sl_k::Ptr{T}
    sl_size::T
end

mutable struct ScsSolution <: AbstractSCSType
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    s::Ptr{Cdouble}
end

mutable struct ScsInfo{T} <: AbstractSCSType
    iter::T
    status::NTuple{128,Cchar}
    lin_sys_solver::NTuple{128,Cchar}
    status_val::T
    scale_updates::T
    pobj::Cdouble
    dobj::Cdouble
    res_pri::Cdouble
    res_dual::Cdouble
    gap::Cdouble
    res_infeas::Cdouble
    res_unbdd_a::Cdouble
    res_unbdd_p::Cdouble
    setup_time::Cdouble
    solve_time::Cdouble
    scale::Cdouble
    comp_slack::Cdouble
    rejected_accel_steps::T
    accepted_accel_steps::T
    lin_sys_time::Cdouble
    cone_time::Cdouble
    accel_time::Cdouble
    ScsInfo{T}() where {T} = new{T}()
end

mutable struct Solution{T}
    x::Array{Float64,1}
    y::Array{Float64,1}
    s::Array{Float64,1}
    info::ScsInfo{T}
    ret_val::T
end

"""
    _ScsDataWrapper

A type for wrapping all data inputs to SCS to prevent them from being freed by
the garbage collector during a solve.

You should not construct this manually. Call `scs_solve` instead.
"""
struct _ScsDataWrapper{S,T}
    linear_solver::S
    m::T
    n::T
    Avalues::Vector{Float64}
    Arowval::Vector{T}
    Acolptr::Vector{T}
    A::ScsMatrix{T}
    Pvalues::Vector{Float64}
    Prowval::Vector{T}
    Pcolptr::Vector{T}
    P::ScsMatrix{T}
    b::Vector{Cdouble}
    c::Vector{Cdouble}
    z::T
    l::T
    bu::Vector{Cdouble}
    bl::Vector{Cdouble}
    q::Vector{T}
    s::Vector{T}
    cs::Vector{T}
    ep::T
    ed::T
    p::Vector{Cdouble}
    d::Vector{T}
    nuc_m::Vector{T}
    nuc_n::Vector{T}
    ell1::Vector{T}
    sl_n::Vector{T}
    sl_k::Vector{T}
    primal::Vector{Cdouble}
    dual::Vector{Cdouble}
    slack::Vector{Cdouble}
    settings::ScsSettings{T}
    options::Any
end

function _sanitize_options(options)
    option_dict = Dict{Symbol,Any}()
    for (key, value) in options
        if key == :linear_solver
            continue
        elseif !hasfield(ScsSettings, key)
            throw(
                ArgumentError(
                    "Unrecognized option passed to SCS solver: $(key)",
                ),
            )
        end
        option_dict[key] = value
    end
    return option_dict
end

function raw_status(info::ScsInfo)
    data = UInt8[info.status[i] for i in 1:findfirst(iszero, info.status)-1]
    return String(data)
end

function _to_sparse(::Type{T}, A::AbstractMatrix) where {T}
    return _to_sparse(T, SparseArrays.sparse(A))
end

function _to_sparse(::Type{T}, A::SparseArrays.SparseMatrixCSC) where {T}
    rowval = convert(Vector{T}, A.rowval .- T(1))
    colptr = convert(Vector{T}, A.colptr .- T(1))
    return A.nzval, rowval, colptr
end

"""
    scs_solve(linear_solver, args...)

The low-level interface to the SCS solver.

!!! warning
    This function is an advanced feature with a risk of incorrect usage. If you
    are a new user, we recommend that you use the JuMP or Convex interfaces
    instead.

## Problem definition

SCS solves a problem of the form
```
minimize        1/2 * x' * P * x + c' * x
subject to      A * x + s = b
                s in K
```
where K is a product cone of

- zero cone `{ 0 }`
- positive orthant `{ x | x ≥ 0 }`
- box cone `{ (t,x) | t*l ≤ x ≤ t*u }`
- second-order cone (SOC) `{ (t,x) | ||x||_2 ≤ t }`
- positive semi-definite (PSD) cone `{ X ∈ ℝⁿˣⁿ | X is psd }`
- complex positive semi-definite cone `{ X ∈ ℂⁿˣⁿ | X is psd }`
- exponential cone `{ (x,y,z) | y e^(x/y) ≤ z, y > 0 }`
- power cone `{ (x,y,z) | x^a * y^(1-a) ≥ |z|, x ≥ 0, y ≥ 0 }`
- dual power cone `{ (u,v,w) | (u/a)^a * (v/(1-a))^(1-a) ≥ |w|, u ≥ 0, v ≥ 0 }`
- log determinant cone `{ (t, u, X) | t ≥ -u log(det(X / u)) }`
- nuclear norm cone `{ (t, X) | t ≥ ||X||_* }`
- L1-matrix norm cone `{ (t, X) | t ≥ Σᵢ|Xᵢⱼ| ∀j}`
- Ky-Fan norm cone `{ (t, X) | t ≥ Σᵏ λᵢ(X) }`

## Input arguments

The problem data are:

- `linear_solver`: a `LinearSolver` to use
- `m`: the number of affine constraints
- `n`: the number of variables
- `A`: an `AbstractMatrix` with `m` rows and `n` cols
- `P`: the upper-triangular component of a positive semidefinite symmetric
       matrix of size `(n, n)`
- `b`: a `Vector` of length `m`
- `c`: a `Vector` of length `n`

The rows of `A` correspond to cones in `K` and need to be specified in the order
above by the following arguments.

- `z`: the number of primal zero / dual free cones, i.e. primal equality
  constraints
- `l`: the number of linear cones
- `bu`: the `Vector` of upper bounds for the box cone
- `bl`: the `Vector` of lower bounds for the box cone
- `q`: the `Vector` of SOCs sizes
- `s`: the `Vector` of PSD cones sizes
- `cs`: the `Vector` of complex PSD cones sizes
- `ep`: the number of primal exponential cones
- `ed`: the number of dual exponential cones
- `p`: the `Vector` of power cone parameters (±1, with negative values for the
  dual cone)
- `d`: the `Vector` of log-det cone sizes
- `nuc_m`: the `Vector` of nuclear-norm row dimensions
- `nuc_n`: the `Vector` of nuclear-norm column dimensions
- `ell1`: the `Vector` of L1-matrix-norm sizes
- `sl_n`: the `Vector` of Ky-Fan norm cone sizes
- `sl_k`: the `Vector` of Ky-Fan norm cone constants

Provide a warm start to SCS by overriding:
- `primal_sol = zeros(n)`: a `Vector` to warmstart the primal variables,
- `dual_sol = zeros(m)`: a `Vector` to warmstart the dual variables,
- `slack = zeros(m)`: a `Vector` to warmstart the slack variables.

!!! note
    To successfully warmstart the solver `primal_sol`, `dual_sol` and `slack`
    must all be provided **and** `warm_start` option must be set to `true`.

Finally, a number of kwarg arguments may be passed as options to the solver.
The valid keywords are `fieldnames(ScsSettings)`. See the official SCS
documentation for information on each keyword.

## Output

This function eturns a `Solution` object, which contains the following fields:
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

!!! warning
    SCS expects the PSD cones to be scaled by a factor of √2. That is,
    the off-diagonal elements in the A matrix and b vector should be multiplied
    by √2, and then the corresponding rows in the dual and slack solution
    vectors should be multiplied by 1/√2.

`SCS.raw_status(::ScsInfo)::String` describes the status, e.g. `"solved"`,
`"infeasible"`, `"unbounded"`, etc. For the precise return status, the value of
`ret_val` field should be compared with
https://github.com/cvxgrp/scs/blob/3aaa93c7aa04c7001df5e51b81f21b126dfa99b3/include/glbopts.h#L18.
"""
function scs_solve(
    linear_solver::Type{<:LinearSolver},
    m::Integer,
    n::Integer,
    A,
    P,
    b::Vector{Float64},
    c::Vector{Float64},
    z::Integer,
    l::Integer,
    bu::Vector{Float64},
    bl::Vector{Float64},
    q::Vector{<:Integer},
    s::Vector{<:Integer},
    cs::Vector{<:Integer},
    ep::Integer,
    ed::Integer,
    p::Vector{Float64},
    d::Vector{<:Integer},
    nuc_m::Vector{<:Integer},
    nuc_n::Vector{<:Integer},
    ell1::Vector{<:Integer},
    sl_n::Vector{<:Integer},
    sl_k::Vector{<:Integer},
    primal_sol::Vector{Float64} = zeros(n),
    dual_sol::Vector{Float64} = zeros(m),
    slack::Vector{Float64} = zeros(m);
    warm_start::Bool = false,
    options...,
)
    if length(primal_sol) == n && length(dual_sol) == length(slack) == m
        if !warm_start
            fill!(primal_sol, 0.0)
            fill!(dual_sol, 0.0)
            fill!(slack, 0.0)
        end
    else
        if warm_start
            throw(ArgumentError("Warmstart doesn't match the problem size"))
        end
        primal_sol, dual_sol, slack = zeros(n), zeros(m), zeros(m)
    end
    @assert length(sl_n) == length(sl_k)
    @assert length(nuc_m) == length(nuc_n)
    T = scsint_t(linear_solver)
    m, n, z, l, ep, ed = T(m), T(n), T(z), T(l), T(ep), T(ed)
    Avalues, Arowval, Acolptr = _to_sparse(T, A)
    Pvalues, Prowval, Pcolptr = _to_sparse(T, P)
    option_dict = _sanitize_options(options)
    option_dict[:warm_start] = convert(Int, warm_start)
    model = _ScsDataWrapper(
        linear_solver,
        m,
        n,
        Avalues,
        Arowval,
        Acolptr,
        ScsMatrix(pointer(Avalues), pointer(Arowval), pointer(Acolptr), m, n),
        Pvalues,
        Prowval,
        Pcolptr,
        ScsMatrix(pointer(Pvalues), pointer(Prowval), pointer(Pcolptr), n, n),
        b,
        c,
        z,
        l,
        bu,
        bl,
        convert(Vector{T}, q),
        convert(Vector{T}, s),
        convert(Vector{T}, cs),
        ep,
        ed,
        p,
        convert(Vector{T}, d),
        convert(Vector{T}, nuc_m),
        convert(Vector{T}, nuc_n),
        convert(Vector{T}, ell1),
        convert(Vector{T}, sl_n),
        convert(Vector{T}, sl_k),
        primal_sol,
        dual_sol,
        slack,
        ScsSettings(linear_solver),
        option_dict,
    )
    Base.GC.@preserve model begin
        return _unsafe_scs_solve(model)
    end
end

function _unsafe_scs_solve(model::_ScsDataWrapper{S,T}) where {S,T}
    scs_cone = ScsCone{T}(
        model.z,
        model.l,
        pointer(model.bu),
        pointer(model.bl),
        max(length(model.bu) - 1, 0),
        pointer(model.q),
        length(model.q),
        pointer(model.s),
        length(model.s),
        pointer(model.cs),
        length(model.cs),
        model.ep,
        model.ed,
        pointer(model.p),
        length(model.p),
        pointer(model.d),
        length(model.d),
        pointer(model.nuc_m),
        pointer(model.nuc_n),
        length(model.nuc_m),
        pointer(model.ell1),
        length(model.ell1),
        pointer(model.sl_k),
        pointer(model.sl_n),
        length(model.sl_n),
    )
    scs_data = ScsData{T}(
        model.m,
        model.n,
        pointer_from_objref(model.A),
        iszero(model.Pvalues) ? C_NULL : pointer_from_objref(model.P),
        pointer(model.b),
        pointer(model.c),
    )
    for (key, value) in model.options
        if value isa AbstractString
            if value isa String
                cstr =
                    Base.unsafe_convert(Cstring, Base.cconvert(Cstring, value))
                setproperty!(model.settings, key, cstr)
            else
                error("You must pass a `String` as the value for $(key)")
            end
        else
            setproperty!(model.settings, key, value)
        end
    end
    scs_work = scs_init(model.linear_solver, scs_data, scs_cone, model.settings)
    scs_solution = ScsSolution(
        pointer(model.primal),
        pointer(model.dual),
        pointer(model.slack),
    )
    scs_info = ScsInfo{T}()
    status = scs_solve(
        model.linear_solver,
        scs_work,
        scs_solution,
        scs_info,
        model.options[:warm_start],
    )
    scs_finish(model.linear_solver, scs_work)
    # SCS does not flush stdout on exit, and the C stdout is not the same as
    # Julia's stdout, so regular calls to flush(stdout) do not work. A
    # work-around is for us to manually flush the output after each solve.
    Base.Libc.flush_cstdio()
    return Solution(model.primal, model.dual, model.slack, scs_info, status)
end

# The following is a backward compatible method that does not have the `cs`
# argument that was introduced in SCS.jl v2.2.0. The underlying C API had a
# breaking change, but we don't need to expose it to the users of SCS.jl.
#
# It is important the function signature of this method matches the function
# signature of the full `scs_solve` exactly to avoid ambiguity errors.
function scs_solve(
    linear_solver::Type{<:LinearSolver},
    m::Integer,
    n::Integer,
    A,
    P,
    b::Vector{Float64},
    c::Vector{Float64},
    z::Integer,
    l::Integer,
    bu::Vector{Float64},
    bl::Vector{Float64},
    q::Vector{<:Integer},
    s::Vector{<:Integer},
    # cs::Vector{<:Integer},  # Skip this argument
    ep::Integer,
    ed::Integer,
    p::Vector{Float64},
    # d::Vector{<:Integer},     # Skip this argument
    # nuc_m::Vector{<:Integer}, # Skip this argument
    # nuc_n::Vector{<:Integer}, # Skip this argument
    # ell1::Vector{<:Integer},  # Skip this argument
    # sl_n::Vector{<:Integer},  # Skip this argument
    # sl_k::Vector{<:Integer},  # Skip this argument
    primal_sol::Vector{Float64} = zeros(n),
    dual_sol::Vector{Float64} = zeros(m),
    slack::Vector{Float64} = zeros(m);
    warm_start::Bool = false,
    options...,
)
    T = scsint_t(linear_solver)
    return scs_solve(
        linear_solver,
        m,
        n,
        A,
        P,
        b,
        c,
        z,
        l,
        bu,
        bl,
        q,
        s,
        T[],  # Default cs argument
        ep,
        ed,
        p,
        T[],  # Default d argument
        T[],  # Default nuc_m argument
        T[],  # Default nuc_n argument
        T[],  # Default ell1 argument
        T[],  # Default sl_n argument
        T[],  # Default sl_k argument
        primal_sol,
        dual_sol,
        slack;
        warm_start,
        options...,
    )
end

function scs_solve(
    linear_solver::Type{<:LinearSolver},
    m::Integer,
    n::Integer,
    A,
    P,
    b::Vector{Float64},
    c::Vector{Float64},
    z::Integer,
    l::Integer,
    bu::Vector{Float64},
    bl::Vector{Float64},
    q::Vector{<:Integer},
    s::Vector{<:Integer},
    cs::Vector{<:Integer},
    ep::Integer,
    ed::Integer,
    p::Vector{Float64},
    # d::Vector{<:Integer},     # Skip this argument
    # nuc_m::Vector{<:Integer}, # Skip this argument
    # nuc_n::Vector{<:Integer}, # Skip this argument
    # ell1::Vector{<:Integer},  # Skip this argument
    # sl_n::Vector{<:Integer},  # Skip this argument
    # sl_k::Vector{<:Integer},  # Skip this argument
    primal_sol::Vector{Float64} = zeros(n),
    dual_sol::Vector{Float64} = zeros(m),
    slack::Vector{Float64} = zeros(m);
    warm_start::Bool = false,
    options...,
)
    T = scsint_t(linear_solver)
    return scs_solve(
        linear_solver,
        m,
        n,
        A,
        P,
        b,
        c,
        z,
        l,
        bu,
        bl,
        q,
        s,
        cs,
        ep,
        ed,
        p,
        T[],  # Default d argument
        T[],  # Default nuc_m argument
        T[],  # Default nuc_n argument
        T[],  # Default ell1 argument
        T[],  # Default sl_n argument
        T[],  # Default sl_k argument
        primal_sol,
        dual_sol,
        slack;
        warm_start,
        options...,
    )
end
