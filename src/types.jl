using SparseArrays

export SCSMatrix, SCSData, SCSSettings, SCSSolution, SCSInfo, SCSCone

VecOrMatOrSparse = Union{VecOrMat, SparseMatrixCSC{Float64,Int}}

SCSInt = Union{Int32, Int64}
abstract type LinearSolver end
struct DirectSolver <: LinearSolver end
struct IndirectSolver <: LinearSolver end
struct IndirectGpuSolver <: LinearSolver end

scsint_t(::Type{<:LinearSolver}) = Int
scsint_t(::Type{IndirectGpuSolver}) = Int32

clib(::Type{DirectSolver}) = direct
clib(::Type{IndirectSolver}) = indirect
clib(::Type{IndirectGpuSolver}) = indirectgpu


struct SCSMatrix{T<:SCSInt}
    values::Ptr{Cdouble}
    rowval::Ptr{T}
    colptr::Ptr{T}
    m::T
    n::T
end

# Version where Julia manages the memory for the vectors.
struct ManagedSCSMatrix{T<:SCSInt}
    values::Vector{Cdouble}
    rowval::Vector{T}
    colptr::Vector{T}
    scsmatref::Base.RefValue{SCSMatrix{T}}

    function ManagedSCSMatrix{T}(m::Integer, n::Integer, A::VecOrMatOrSparse) where T
        A_sparse = sparse(A)

        # we convert everything to make sure that no conversion will
        # accidentally happen when ccalling `new`, so that `scsmatref`
        # holds pointers to actual data stored in the fields.

        values = convert(Vector{Cdouble}, copy(A_sparse.nzval))
        rowval = convert(Vector{T}, A_sparse.rowval .- 1)
        colptr = convert(Vector{T}, A_sparse.colptr .- 1)

        # scsmatref holds the reference to SCSMatrix created out of data in ManagedSCSMatrix.
        # this way the reference to SCSMatrix always points to valid data
        # as long as the ManagedSCSMatrix is not GC collected.
        # One MUST
        #   `Base.unsafe_convert(Ref{SCSMatrix}, scsmatref)`
        # when assigning to a field of type `Ptr{SCSMatrix}`, or specify
        #   `Ref{SCSMatrix}`
        # in the type tuple when ccalling with `scsmatref`

        scsmat = SCSMatrix{T}(
            pointer(values), pointer(rowval), pointer(colptr), m, n)

        return new{T}(values, rowval, colptr, Base.cconvert(Ref{SCSMatrix{T}}, scsmat))
    end
end

struct SCSSettings{T<:SCSInt}
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

    SCSSettings{T}() where T = new{T}()
    SCSSettings{T}(normalize, scale, rho_x, max_iters, eps, alpha, cg_rate, verbose, warm_start, acceleration_lookback, write_data_filename) where T =
        new{T}(normalize, scale, rho_x, max_iters, eps, alpha, cg_rate, verbose, warm_start, acceleration_lookback, write_data_filename)
end

function _SCS_user_settings(default_settings::SCSSettings{T};
        normalize=default_settings.normalize,
        scale=default_settings.scale,
        rho_x=default_settings.rho_x,
        max_iters=default_settings.max_iters,
        eps=default_settings.eps,
        alpha=default_settings.alpha,
        cg_rate=default_settings.cg_rate,
        verbose=default_settings.verbose,
        warm_start=default_settings.warm_start,
        acceleration_lookback=default_settings.acceleration_lookback,
        write_data_filename=default_settings.write_data_filename
        ) where T
    return SCSSettings{T}(normalize, scale, rho_x, max_iters, eps, alpha, cg_rate, verbose,warm_start, acceleration_lookback, write_data_filename)
end

function SCSSettings(linear_solver::Type{<:LinearSolver}; options...)

    T = scsint_t(linear_solver)

    managed_matrix = ManagedSCSMatrix{T}(0,0,spzeros(0,0))
    default_settings = Base.cconvert(Ref{SCSSettings{T}}, SCSSettings{T}())
    a = [0.0]
    dummy_data = SCSData(0,0,
        managed_matrix,
        a,
        a,
        default_settings)

    Base.GC.@preserve managed_matrix default_settings a begin
        SCS_set_default_settings(linear_solver, dummy_data)
    end

    return _SCS_user_settings(default_settings[]; options...)
end

struct SCSData{T<:SCSInt}
    # A has m rows, n cols
    m::T
    n::T
    A::Ptr{SCSMatrix{T}}
    # b is of size m, c is of size n
    b::Ptr{Cdouble}
    c::Ptr{Cdouble}
    stgs::Ptr{SCSSettings{T}}
end

function SCSData(m::Integer, n::Integer,
    mat::ManagedSCSMatrix{T},
    b::Vector{Float64},
    c::Vector{Float64},
    stgs::Ref{SCSSettings{T}}) where T
    return SCSData{T}(m, n,
        Base.unsafe_convert(Ref{SCSMatrix{T}}, mat.scsmatref), # Ptr{SCSMatrix{T}}
        pointer(b),
        pointer(c),
        Base.unsafe_convert(Ref{SCSSettings{T}}, stgs) # Ptr{SCSSettings{T}}
    )
end

struct SCSSolution
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    s::Ptr{Cdouble}
end

struct SCSInfo{T<:SCSInt}
    iter::T
    status::NTuple{32, Cchar} # char status[32]
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
end

SCSInfo{T}() where T = SCSInfo{T}(0, ntuple(_ -> zero(Cchar), 32), 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

function raw_status(info::SCSInfo)
    s = collect(info.status)
    len = findfirst(iszero, s) - 1
    # There is no method String(::Vector{Cchar}) so we convert to `UInt8`.
    return String(UInt8[s[i] for i in 1:len])
end

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
struct SCSCone{T<:SCSInt}
    f::T # number of linear equality constraints
    l::T # length of LP cone
    q::Ptr{T} # array of second-order cone constraints
    qsize::T # length of SOC array
    s::Ptr{T} # array of SD constraints
    ssize::T # length of SD array
    ep::T # number of primal exponential cone triples
    ed::T # number of dual exponential cone triples
    p::Ptr{Cdouble} # array of power cone params, must be \in [-1, 1], negative values are interpreted as specifying the dual cone
    psize::T # length of power cone array
end

# Returns an SCSCone. The q, s, and p arrays are *not* GC tracked in the
# struct. Use this only when you know that q, s, and p will outlive the struct.
function SCSCone{T}(f::Integer, l::Integer, q::Vector{T}, s::Vector{T},
                 ep::Integer, ed::Integer, p::Vector{Cdouble}) where T
    return SCSCone{T}(f, l, pointer(q), length(q), pointer(s), length(s), ep, ed, pointer(p), length(p))
end


mutable struct Solution{T<:SCSInt}
    x::Array{Float64, 1}
    y::Array{Float64, 1}
    s::Array{Float64, 1}
    info::SCSInfo{T}
    ret_val::T
end

function sanitize_SCS_options(options)
    options = Dict(options)
    if haskey(options, :linear_solver)
        linear_solver = options[:linear_solver]
        if !(linear_solver in available_solvers)
            throw(ArgumentError("Unrecognized linear_solver passed to SCS: $linear_solver;\nRecognized options are: $(join(available_solvers, ", ", " and "))."))
        end
        delete!(options, :linear_solver)
    else
        linear_solver = IndirectSolver # the default linear_solver
    end

    SCS_options = append!([:linear_solver], fieldnames(SCSSettings))
    unrecognized = setdiff(keys(options), SCS_options)
    if length(unrecognized) > 0
        plur = length(unrecognized) > 1 ? "s" : ""
        throw(ArgumentError("Unrecognized option$plur passed to SCS: $(join(unrecognized, ", "));\nRecognized options are: $(join(SCS_options, ", ", " and "))."))
    end
    return linear_solver, options
end
