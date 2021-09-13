abstract type LinearSolver end

abstract type AbstractSCSType end

Base.cconvert(::Type{Ptr{Cvoid}}, x::AbstractSCSType) = x

function Base.unsafe_convert(::Type{Ptr{Cvoid}}, x::AbstractSCSType)
    return pointer_from_objref(x)
end

mutable struct SCSMatrix{T} <: AbstractSCSType
    values::Ptr{Cdouble}
    rowval::Ptr{T}
    colptr::Ptr{T}
    m::T
    n::T
end

mutable struct SCSSettings{T} <: AbstractSCSType
    normalize::T # boolean, heuristic data rescaling
    scale::Cdouble # if normalized, rescales by this factor
    rho_x::Cdouble # x equality constraint scaling
    max_iters::T # maximum iterations to take
    eps::Cdouble # convergence tolerance
    alpha::Cdouble # relaxation parameter
    cg_rate::Cdouble # for indirect, tolerance goes down like (1/iter)^cg_rate
    verbose::T # boolean, write out progress
    warm_start::T # boolean, warm start (put initial guess in Sol struct)
    acceleration_lookback::T # acceleration memory parameter
    write_data_filename::Cstring
    SCSSettings{T}() where {T} = new{T}()
end

mutable struct SCSData{T} <: AbstractSCSType
    m::T
    n::T
    A::Ptr{Cvoid}
    b::Ptr{Cdouble}
    c::Ptr{Cdouble}
    stgs::Ptr{Cvoid}
end

mutable struct SCSSolution <: AbstractSCSType
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    s::Ptr{Cdouble}
end

mutable struct SCSInfo{T} <: AbstractSCSType
    iter::T
    status::NTuple{32,Cchar}
    statusVal::T
    pobj::Cdouble
    dobj::Cdouble
    resPri::Cdouble
    resDual::Cdouble
    resInfeas::Cdouble
    resUnbdd::Cdouble
    relGap::Cdouble
    setupTime::Cdouble
    solveTime::Cdouble
    SCSInfo{T}() where {T} = new{T}()
end

mutable struct SCSCone{T} <: AbstractSCSType
    f::T
    l::T
    q::Ptr{T}
    qsize::T
    s::Ptr{T}
    ssize::T
    ep::T
    ed::T
    p::Ptr{Cdouble}
    psize::T
end

mutable struct Solution{T}
    x::Array{Float64,1}
    y::Array{Float64,1}
    s::Array{Float64,1}
    info::SCSInfo{T}
    ret_val::T
end

"""
    _SCSDataWrapper

A type for wrapping all data inputs to SCS to prevent them from being freed by
the garbage collector during a solve.

You should not construct this manually. Call `SCS_solve` instead.
"""
struct _SCSDataWrapper{S,T}
    linear_solver::S
    m::T
    n::T
    A::SCSMatrix{T}
    values::Vector{Cdouble}
    rowval::Vector{T}
    colptr::Vector{T}
    b::Vector{Cdouble}
    c::Vector{Cdouble}
    f::T
    l::T
    q::Vector{T}
    s::Vector{T}
    ep::T
    ed::T
    p::Vector{Cdouble}
    primal::Vector{Cdouble}
    dual::Vector{Cdouble}
    slack::Vector{Cdouble}
    settings::SCSSettings
    options::Any
end

function _sanitize_options(options)
    option_dict = Dict{Symbol,Any}()
    fields = fieldnames(SCSSettings)
    for (key, value) in options
        if key == :linear_solver
            continue
        elseif !(key in fields)
            throw(ArgumentError("Unrecognized option $(key)"))
        end
        option_dict[key] = value
    end
    return option_dict
end

function raw_status(info::SCSInfo)
    data = UInt8[info.status[i] for i in 1:findfirst(iszero, info.status)-1]
    return String(data)
end

function _to_sparse(::Type{T}, A::AbstractMatrix) where {T}
    sparse_A = SparseArrays.sparse(A)
    values = sparse_A.nzval
    rowval = convert(Vector{T}, sparse_A.rowval .- 1)
    colptr = convert(Vector{T}, sparse_A.colptr .- 1)
    return values, rowval, colptr
end

"""
    SCS_solve(args...)

SCS solves a problem of the form
```
minimize        c' * x
subject to      A * x + s = b
                s in K
```
where K is a product cone of
 * zero cone { x | x = 0 },
 * positive orthant { x | x >= 0 },
 * second-order cones (SOC) { (t,x) | ||x||_2 <= t },
 * semi-definite cones (SDC) { X | X psd }, and
 * exponential cones {(x,y,z) | y e^(x/y) <= z, y>0 },
 * power cone { (x,y,z) | x^a*y^(1-a) >= |z|, x>=0, y>=0 },
 * dual power cone { (u,v,w) | (u/a)^a * (v/(1-a))^(1-a) >= |w|, u>=0, v>=0}

Description of input argments:
 * A is the matrix with m rows and n cols
 * b is a vector of length m
 * c is a vector of length n

The rows of A correspond to cones in K and need to be specified in the order
above.

 * f (num primal zero / dual free cones, i.e. primal equality constraints)
 * l (num linear cones)
 * q (array of SOCs sizes)
 * s (array of SDCs sizes)
 * ep (num primal exponential cones)
 * ed (num dual exponential cones)
 * p (array of power cone params, must be in [-1, 1], negative values specify
   the dual cone)

Returns a Solution object.
"""
function SCS_solve(
    linear_solver::Type{<:LinearSolver},
    m::Integer,
    n::Integer,
    A::AbstractMatrix,
    b::Vector{Float64},
    c::Vector{Float64},
    f::Integer,
    l::Integer,
    q::Vector{<:Integer},
    s::Vector{<:Integer},
    ep::Integer,
    ed::Integer,
    p::Vector{Float64},
    primal_sol::Vector{Float64} = zeros(n),
    dual_sol::Vector{Float64} = zeros(m),
    slack::Vector{Float64} = zeros(m);
    warm_start::Bool = false,
    options...,
)
    if n <= 0
        throw(ArgumentError("The number of variables must be greater than 0"))
    elseif m <= 0
        throw(ArgumentError("The number of constraints must be greater than 0"))
    end
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
    T = scsint_t(linear_solver)
    m, n, f, l, ep, ed = T(m), T(n), T(f), T(l), T(ep), T(ed)
    values, rowval, colptr = _to_sparse(T, A)
    option_dict = _sanitize_options(options)
    if warm_start
        option_dict[:warm_start] = 1
    end
    model = _SCSDataWrapper(
        linear_solver,
        m,
        n,
        SCSMatrix(pointer(values), pointer(rowval), pointer(colptr), m, n),
        values,
        rowval,
        colptr,
        b,
        c,
        f,
        l,
        convert(Vector{T}, q),
        convert(Vector{T}, s),
        ep,
        ed,
        p,
        primal_sol,
        dual_sol,
        slack,
        SCSSettings{T}(),
        option_dict,
    )
    Base.GC.@preserve model begin
        return _unsafe_scs_solve(model)
    end
end

function _unsafe_scs_solve(model::_SCSDataWrapper{S,T}) where {S,T}
    scs_solution = SCSSolution(
        pointer(model.primal),
        pointer(model.dual),
        pointer(model.slack),
    )
    scs_cone = SCSCone{T}(
        model.f,
        model.l,
        pointer(model.q),
        length(model.q),
        pointer(model.s),
        length(model.s),
        model.ep,
        model.ed,
        pointer(model.p),
        length(model.p),
    )
    scs_info = SCSInfo{T}()
    scs_data = SCSData{T}(
        model.m,
        model.n,
        pointer_from_objref(model.A),
        pointer(model.b),
        pointer(model.c),
        pointer_from_objref(model.settings),
    )
    SCS_set_default_settings(model.linear_solver, scs_data)
    for (key, value) in model.options
        setfield!(model.settings, key, value)
    end
    p = SCS_init(model.linear_solver, scs_data, scs_cone, scs_info)
    status = SCS_solve(
        model.linear_solver,
        p,
        scs_data,
        scs_cone,
        scs_solution,
        scs_info,
    )
    SCS_finish(model.linear_solver, p)
    return Solution(model.primal, model.dual, model.slack, scs_info, status)
end
