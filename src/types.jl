using Compat
export SCSMatrix, SCSData, SCSSolution, SCSInfo, SCSCone, SCSWork, SCSVecOrMatOrSparse


SCSVecOrMatOrSparse = Union(VecOrMat, SparseMatrixCSC{Float64,Int})


immutable SCSMatrix
    m::Int
    n::Int
    values::Ptr{Cdouble}
    rowval::Ptr{Int}
    colptr::Ptr{Int}
end


immutable SCSData
    # A has m rows, n cols
    m::Int
    n::Int
    A::Ptr{SCSMatrix}
    # b is of size m, c is of size n
    b::Ptr{Cdouble}
    c::Ptr{Cdouble}
    # max_iters to take
    max_iters::Int
    # convergence tolerance
    eps::Cdouble
    # relaxation parameter
    alpha::Cdouble
    # x equality constraint scaling
    rho_x::Cdouble
    # for indirect, tolerance goes down like (1/iter)^CG_RATE: 1.5
    cg_rate::Cdouble
    verbose::Int
    # 0 or 1
    normalize::Int
    # if normalized, rescales by this factor
    scale::Cdouble
    # 0 or 1
    warm_start::Int
end


immutable SCSSolution
    x::Ptr{None}
    y::Ptr{None}
    s::Ptr{None}
end


immutable SCSInfo
    iter::Int
    # We need to allocate 32 bytes for a character string, so we allocate 256 bits
    # of integer instead
    # TODO: Find a better way to do this
    status1::Int128
    status2::Int128

    statusVal::Int
    pobj::Cdouble
    dobj::Cdouble
    resPri::Cdouble
    resDual::Cdouble
    resInfeas::Cdouble
    resUnbdd::Cdouble
    relGap::Cdouble
    setupTime::Cdouble
    solveTime::Cdouble
end


immutable SCSCone
    f::Int # number of linear equality constraints
    l::Int # length of LP cone
    q::Ptr{Int} # array of second-order cone constraints
    qsize::Int # length of SOC array
    s::Ptr{Int} # array of semi-definite constraints
    ssize::Int # length of SD array
    ep::Int # number of primal exponential cone triples
    ed::Int # number of dual exponential cone triples
end


immutable SCSWork
    u::Ptr{Cdouble}
    v::Ptr{Cdouble}
    u_t::Ptr{Cdouble}
    u_prev::Ptr{Cdouble}
    h::Ptr{Cdouble}
    g::Ptr{Cdouble}
    pr::Ptr{Cdouble}
    dr::Ptr{Cdouble}
    gTh::Cdouble
    sc_b::Cdouble
    sc_c::Cdouble
    nm_b::Cdouble
    nm_c::Cdouble
    meanNormRowA::Cdouble
    meanNormColA::Cdouble
    D::Ptr{Cdouble}
    E::Ptr{Cdouble}
    p::Ptr{Void}
end

@compat const status_map = Dict{Int, Symbol}(
    1 => :Optimal,
    -2 => :Infeasible,
    -1 => :Unbounded,
    -3 => :Indeterminate,
    -4 => :Error
)

type Solution
    x::Array{Float64, 1}
    y::Array{Float64, 1}
    s::Array{Float64, 1}
    status::Symbol
    ret_val::Int

    function Solution(x::Array{Float64, 1}, y::Array{Float64, 1}, s::Array{Float64, 1}, ret_val::Int)
        if haskey(status_map, ret_val)
            return new(x, y, s, status_map[ret_val], ret_val)
        else
            return new(x, y, s, :UnknownError, ret_val)
        end
    end
end
