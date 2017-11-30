export SCSMatrix, SCSData, SCSSettings, SCSSolution, SCSInfo, SCSCone, SCSVecOrMatOrSparse


SCSVecOrMatOrSparse = Union{VecOrMat, SparseMatrixCSC{Float64,Int}}


immutable SCSMatrix
    values::Ptr{Cdouble}
    rowval::Ptr{Int}
    colptr::Ptr{Int}
    m::Int
    n::Int
end


immutable SCSSettings
    normalize::Int # boolean, heuristic data rescaling: 1
    scale::Cdouble # if normalized, rescales by this factor: 1
    rho_x::Cdouble # x equality constraint scaling: 1e-3
    max_iters::Int # maximum iterations to take: 2500
    eps::Cdouble # convergence tolerance: 1e-3
    alpha::Cdouble # relaxation parameter: 1.8
    cg_rate::Cdouble # for indirect, tolerance goes down like (1/iter)^cg_rate: 2
    verbose::Int # boolean, write out progress: 1
    warm_start::Int # boolean, warm start (put initial guess in Sol struct): 0
    acceleration_lookback::Int # boolean, acceleration memory parameter: 20
end


immutable SCSData
    # A has m rows, n cols
    m::Int
    n::Int
    A::Ptr{SCSMatrix}
    # b is of size m, c is of size n
    b::Ptr{Cdouble}
    c::Ptr{Cdouble}
    stgs::Ptr{SCSSettings}
end


immutable SCSSolution
    x::Ptr{Void}
    y::Ptr{Void}
    s::Ptr{Void}
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
    p::Ptr{Cdouble} # array of power cone params
    psize::Int # length of power cone array
end


# TODO needs to be updated for newest constants
const status_map = Dict{Int, Symbol}(
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
