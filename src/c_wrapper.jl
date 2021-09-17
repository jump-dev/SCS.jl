abstract type LinearSolver end

abstract type AbstractSCSType end

Base.cconvert(::Type{Ptr{Cvoid}}, x::AbstractSCSType) = x

function Base.unsafe_convert(::Type{Ptr{Cvoid}}, x::AbstractSCSType)
    return pointer_from_objref(x)
end

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
        stgs = new{scsint_t(linear_solver)}()
        scs_set_default_settings(linear_solver, stgs)
        return stgs
    end
end

mutable struct ScsData{T} <: AbstractSCSType
    m::T
    n::T
    A::Ptr{Cvoid}
    P::Ptr{Cvoid}
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
    ep::T
    ed::T
    p::Ptr{Cdouble}
    psize::T
end

mutable struct ScsSolution <: AbstractSCSType
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    s::Ptr{Cdouble}
end

mutable struct ScsInfo{T} <: AbstractSCSType
    iter::T
    status::NTuple{128,Cchar}
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
    values_A::Vector{Float64}
    rowval_A::Vector{T}
    colptr_A::Vector{T}
    A::ScsMatrix{T}
    values_P::Vector{Float64}
    rowval_P::Vector{T}
    colptr_P::Vector{T}
    P::ScsMatrix{T}
    b::Vector{Cdouble}
    c::Vector{Cdouble}
    z::T
    l::T
    bu::Vector{Cdouble}
    bl::Vector{Cdouble}
    q::Vector{T}
    s::Vector{T}
    ep::T
    ed::T
    p::Vector{Cdouble}
    primal::Vector{Cdouble}
    dual::Vector{Cdouble}
    slack::Vector{Cdouble}
    settings::ScsSettings
    options::Any
end

function _sanitize_options(options)
    option_dict = Dict{Symbol,Any}()
    fields = fieldnames(ScsSettings)
    for (key, value) in options
        if key == :linear_solver
            continue
        elseif !(key in fields)
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

SCS solves a problem of the form
```
minimize        1/2 * x' * P * x + c' * x
subject to      A * x + s = b
                s in K
```
where K is a product cone of
- zero cone,
- positive orthant `{ x | x ≥ 0 }`,
- box cone `{ (t,x) | t*l ≤ x ≤ t*u}`,
- second-order cone (SOC) `{ (t,x) | ||x||_2 ≤ t }`,
- semi-definite cone (SDC) `{ X | X is psd }`,
- exponential cone `{ (x,y,z) | y e^(x/y) ≤ z, y>0 }`,
- power cone `{ (x,y,z) | x^a * y^(1-a) ≥ |z|, x ≥ 0, y ≥ 0 }`, and
- dual power cone `{ (u,v,w) | (u/a)^a * (v/(1-a))^(1-a) ≥ |w|, u ≥ 0, v ≥ 0 }`.

The problem data are:
- `A`: `AbstractMatrix` with `m` rows and `n` cols
- `P`: positive semidefinite symmetric `AbstractMatrix` of size `(n, n)`
- `b`: `AbstractVector` of length `m`
- `c`: `AbstractVector` of length `n`.

The rows of `A` correspond to cones in `K` and need to be specified in the order
above by the following arguments.

- `z`: the number of primal zero / dual free cones, i.e. primal equality
  constraints
- `l`: the number of linear cones
- `bu`: the `Vector` of upper bounds for the box cone
- `bl`: the `Vector` of lower bounds for the box cone
- `q`: the `Vector` of SOCs sizes
- `s`: the `Vector` of SDCs sizes
- `ep`: the number of primal exponential cones
- `ed`: the number of dual exponential cones
- `p`: the `Vector` of power cone parameters (±1, with negative values for the
  dual cone)

Additionally one may provide warm start to the solver by overriding
- `primal_sol = zeros(n)`: a `Vector` to warmstart the primal variables,
- `dual_sol = zeros(m)`: a `Vector` to warmstart the dual variables,
- `slack = zeros(m)`: a `Vector` to warmstart the slack variables.

!!! note
    To successfully warmstart the solver all `primal_sol`, `dual_sol` and
    `slack` must be provided **and** `warms_start` option must be set to `true`.

Finally, a number of kwarg arguments may be passed as options to the solver.
The valid keywords are `fieldnames(ScsSettings)`. For the meaning please consult
`scs` documentation.

Returns a Solution object.

!!! warning
    SCS expects the semi-definite cones to be scaled by a factor of √2. That is,
    the off-diagonal elements in the A matrix and b vector should be multiplied
    by √2, and then the corresponding rows in the dual and slack solution
    vectors should be multiplied by 1/√2.
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

    # m == size(A, 1) || throw(
    #     ArgumentError(
    #         "The number of rows in A inconsistent with input dimension",
    #     ),
    # )
    # n == size(A, 2) || throw(
    #     ArgumentError(
    #         "The number of cols in A inconsistent with input dimension",
    #     ),
    # )
    # size(P, 1) == size(P, 2) || throw(ArgumentError("P is not square"))
    # size(P, 1) == n || throw(
    #     ArgumentError("P dimension $(size(P, 1)) inconsistent with n = $n"),
    # )

    length(bu) == length(bl) || throw(
        ArgumentError(
            "The upper and lower bounds of the box cone differ in dimension",
        ),
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
    T = scsint_t(linear_solver)
    m, n, z, l, ep, ed = T(m), T(n), T(z), T(l), T(ep), T(ed)
    values_A, rowval_A, colptr_A = _to_sparse(T, A)
    values_P, rowval_P, colptr_P = _to_sparse(T, P)
    option_dict = _sanitize_options(options)
    if warm_start
        option_dict[:warm_start] = 1
    end
    model = _ScsDataWrapper(
        linear_solver,
        m,
        n,
        values_A,
        rowval_A,
        colptr_A,
        ScsMatrix(
            pointer(values_A),
            pointer(rowval_A),
            pointer(colptr_A),
            m,
            n,
        ),
        values_P,
        rowval_P,
        colptr_P,
        ScsMatrix(
            pointer(values_P),
            pointer(rowval_P),
            pointer(colptr_P),
            n,
            n,
        ),
        b,
        c,
        z,
        l,
        bu,
        bl,
        convert(Vector{T}, q),
        convert(Vector{T}, s),
        ep,
        ed,
        p,
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
    scs_solution = ScsSolution(
        pointer(model.primal),
        pointer(model.dual),
        pointer(model.slack),
    )
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
        model.ep,
        model.ed,
        pointer(model.p),
        length(model.p),
    )
    scs_info = ScsInfo{T}()
    scs_data = ScsData{T}(
        model.m,
        model.n,
        pointer_from_objref(model.A),
        iszero(model.values_P) ? convert(Ptr{Float64}, C_NULL) :
        pointer_from_objref(model.P),
        pointer(model.b),
        pointer(model.c),
    )
    for (key, value) in model.options
        setproperty!(model.settings, key, value)
    end

    @assert scs_data.m == model.A.m
    @assert scs_data.n == model.A.n

    scs_work = scs_init(model.linear_solver, scs_data, scs_cone, model.settings)
    status = scs_solve(model.linear_solver, scs_work, scs_solution, scs_info)
    scs_finish(model.linear_solver, scs_work)
    return Solution(model.primal, model.dual, model.slack, scs_info, status)
end
