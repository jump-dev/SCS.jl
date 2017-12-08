export SCSMatrix, SCSData, SCSSettings, SCSSolution, SCSInfo, SCSCone, SCSVecOrMatOrSparse


SCSVecOrMatOrSparse = Union{VecOrMat, SparseMatrixCSC{Float64,Int}}

struct SCSMatrix
    values::Ptr{Cdouble}
    rowval::Ptr{Int}
    colptr::Ptr{Int}
    m::Int
    n::Int
end

# Version where Julia manages the memory for the vectors.
struct ManagedSCSMatrix
    values::Vector{Cdouble}
    rowval::Vector{Int}
    colptr::Vector{Int}
    m::Int
    n::Int
end

function ManagedSCSMatrix(m::Int, n::Int, A::SCSVecOrMatOrSparse)
    A_sparse = sparse(A)

    values = copy(A_sparse.nzval)
    rowval = convert(Array{Int, 1}, A_sparse.rowval .- 1)
    colptr = convert(Array{Int, 1}, A_sparse.colptr .- 1)

    return ManagedSCSMatrix(values, rowval, colptr, m, n)
end

# Returns an SCSMatrix. The vectors are *not* GC tracked in the struct.
# Use this only when you know that the managed matrix will outlive the SCSMatrix.
SCSMatrix(m::ManagedSCSMatrix) =
    SCSMatrix(pointer(m.values), pointer(m.rowval), pointer(m.colptr), m.m, m.n)


struct SCSSettings
    normalize::Int # boolean, heuristic data rescaling: 1
    scale::Cdouble # if normalized, rescales by this factor: 1
    rho_x::Cdouble # x equality constraint scaling: 1e-3
    max_iters::Int # maximum iterations to take: 5000
    eps::Cdouble # convergence tolerance: 1e-5
    alpha::Cdouble # relaxation parameter: 1.5
    cg_rate::Cdouble # for indirect, tolerance goes down like (1/iter)^cg_rate: 2
    verbose::Int # boolean, write out progress: 1
    warm_start::Int # boolean, warm start (put initial guess in Sol struct): 0
    acceleration_lookback::Int # acceleration memory parameter: 20
end

function SCSSettings(;normalize=1::Int, scale=convert(Cdouble, 1.0)::Cdouble, rho_x=convert(Cdouble,1e-3)::Cdouble,
                        max_iters=5000::Int, eps=convert(Cdouble, 1e-5)::Cdouble, alpha=convert(Cdouble, 1.5)::Cdouble,
                        cg_rate=convert(Cdouble,2)::Cdouble, verbose=1::Int, warm_start=0::Int, acceleration_lookback=20::Int)
    return SCSSettings(normalize, scale, rho_x, max_iters, eps, alpha, cg_rate, verbose, warm_start, acceleration_lookback)
end


struct SCSData
    # A has m rows, n cols
    m::Int
    n::Int
    A::Ptr{SCSMatrix}
    # b is of size m, c is of size n
    b::Ptr{Cdouble}
    c::Ptr{Cdouble}
    stgs::Ptr{SCSSettings}
end

struct SCSSolution
    x::Ptr{Void}
    y::Ptr{Void}
    s::Ptr{Void}
end


struct SCSInfo
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

SCSInfo() = SCSInfo(0, convert(Int128, 0), convert(Int128, 0), 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

# SCS solves a problem of the form
# minimize        c' * x
# subject to      A * x + s = b
#                 s in K
# where K is a product cone of
# zero cones,
# linear cones { x | x >= 0 },
# second-order cones { (t,x) | ||x||_2 <= t },
# semi-definite cones { X | X psd }, and
# exponential cones {(x,y,z) | y e^(x/y) <= z, y>0 }.
# dual exponential cones {(u,v,w) | −u e^(v/u) <= e w, u<0}
# power cones {(x,y,z) | x^a * y^(1-a) >= |z|, x>=0, y>=0}
# dual power cones {(u,v,w) | (u/a)^a * (v/(1-a))^(1-a) >= |w|, u>=0, v>=0}
struct SCSCone
    f::Int # number of linear equality constraints
    l::Int # length of LP cone
    q::Ptr{Int} # array of second-order cone constraints
    qsize::Int # length of SOC array
    s::Ptr{Int} # array of SD constraints
    ssize::Int # length of SD array
    ep::Int # number of primal exponential cone triples
    ed::Int # number of dual exponential cone triples
    p::Ptr{Cdouble} # array of power cone params, must be \in [-1, 1], negative values are interpreted as specifying the dual cone
    psize::Int # length of power cone array
end

# Returns an SCSCone. The q, s, and p arrays are *not* GC tracked in the
# struct. Use this only when you know that q, s, and p will outlive the struct.
function SCSCone(f::Int, l::Int, q::Vector{Int}, s::Vector{Int},
                 ep::Int, ed::Int, p::Vector{Cdouble})
    return SCSCone(f, l, pointer(q), length(q), pointer(s), length(s), ep, ed, pointer(p), length(p))
end


# TODO needs to be updated for newest constants
const status_map = Dict{Int, Symbol}(
    1 => :Optimal,
    -2 => :Infeasible,
    -1 => :Unbounded,
    -3 => :Indeterminate,
    -4 => :Error
)

mutable struct Solution
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
