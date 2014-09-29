export create_scs_matrix, create_scs_data, create_scs_cone


# Takes a vector or matrix or sparse matrix A and creates an SCSMatrix
function create_scs_matrix(A::SCSVecOrMatOrSparse)
    A_sparse = sparse(A)

    values = A_sparse.nzval * 1.0
    rowval = convert(Array{Clong, 1}, A_sparse.rowval .- 1)
    colptr = convert(Array{Clong, 1}, A_sparse.colptr .- 1)

    return SCSMatrix(pointer(values), pointer(rowval), pointer(colptr))
end


# Create an SCSData type
# We assume we are solving a problem of the form
# minimize        c' * x
# subject to      A * x + s = b
#                 s in K
# A is the matrix with m rows and n cols
# b is of length m x 1
# c is of length n x 1
# refer to create_scs_cone for K
function create_scs_data(;m::Clong=nothing, n::Clong=nothing, A::Ptr{SCSMatrix}=nothing,
        b::Ptr{Cdouble}=nothing,  c::Ptr{Cdouble}=nothing, max_iters=2500::Clong,
        eps=convert(Cdouble, 1e-5)::Cdouble, alpha=convert(Cdouble, 1.8)::Cdouble,
        rho_x=convert(Cdouble, 1e-3)::Cdouble, scale=convert(Cdouble, 5.0)::Cdouble,
        cg_rate=convert(Cdouble, 1.5)::Cdouble, verbose=1::Clong,
        normalize=1::Clong, warm_start=0::Clong)

    data = SCSData(m, n, A, b, c, max_iters, eps, alpha, rho_x, cg_rate, verbose, normalize, scale, warm_start)
    return data
end


# Refer to comment above
function create_scs_data(m::Clong, n::Clong, A::Ptr{SCSMatrix}, b::Ptr{Cdouble}, c::Ptr{Cdouble})
    return create_scs_data(m=m, n=n, A=A, b=b, c=c)
end


# Refer to comment above
function create_scs_data(m::Int64, n::Int64, A::SCSVecOrMatOrSparse, b::Array{Float64,},
        c::Array{Float64,})
    if size(b, 1) != m || size(b, 2) != 1 || size(c, 1) != n || size(c, 2) != 1
        error("Size of b must be m x 1 and size of c must be n x 1")
    end
    A = [create_scs_matrix(A)]
    return create_scs_data(m=m, n=n, A=pointer(A), b=pointer(b), c=pointer(c))
end


# Create an SCSCone type
# We assume we are solving a problem of the form
# minimize        c' * x
# subject to      A * x + s = b
#                 s in K
# where K is a product cone of
# zero cones,
# linear cones { x | x >= 0 },
# second-order cones { (t,x) | ||x||_2 <= t },
# semi-definite cones { X | X psd }, and
# exponential cones {(x,y,z) | y e^(x/y) <= z, y>0 }.
#
# f (num primal zero / dual free cones, i.e. primal equality constraints)
# l (num linear cones)
# q (array of SOCs sizes)
# s (array of SDCs sizes)
# ep (num primal exponential cones)
# ed (num dual exponential cones).
function create_scs_cone(f::Clong, l::Clong, q::Ptr{Clong}, qsize::Clong, s::Ptr{Clong},
        ssize::Clong, ep::Clong, ed::Clong)
    return SCSCone(f, l, q, qsize, s, ssize, ep, ed)
end


# Refer to comment above
function create_scs_cone(f::Clong, l::Clong, q::Array{Int64,}, qsize::Clong, s::Array{Int64,},
        ssize::Clong, ep::Clong, ed::Clong)
    return SCSCone(f, l, pointer(q), qsize, pointer(s), ssize, ep, ed)
end


# We assume we are solving a problem of the form
# minimize        c' * x
# subject to      A * x + s = b
#                 s in K
# where K is a product cone of
# zero cones,
# linear cones { x | x >= 0 },
# second-order cones { (t,x) | ||x||_2 <= t },
# semi-definite cones { X | X psd }, and
# exponential cones {(x,y,z) | y e^(x/y) <= z, y>0 }.
#
#
# A is the matrix with m rows and n cols
# b is of length m x 1
# c is of length n x 1
#
# f (num primal zero / dual free cones, i.e. primal equality constraints)
# l (num linear cones)
# q (array of SOCs sizes)
# s (array of SDCs sizes)
# ep (num primal exponential cones)
# ed (num dual exponential cones).
#
# Returns object of type Solution
# type Solution with
# x, y, s, status (ASCII string), ret_val (numerical status)
#
function SCS_solve(m::Int64, n::Int64, A::SCSVecOrMatOrSparse, b::Array{Float64,},
        c::Array{Float64,}, f::Clong, l::Clong, q::Array{Int64,}, qsize::Clong, s::Array{Int64,},
        ssize::Clong, ep::Clong, ed::Clong)
    data = create_scs_data(m, n, A, b, c)
    cone = create_scs_cone(f, l, q, qsize, s, ssize, ep, ed)

    status, solution, info, p_work = SCS_solve(data, cone)

    ptr_x = convert(Ptr{Float64}, solution.x)
    ptr_y = convert(Ptr{Float64}, solution.y)
    ptr_s = convert(Ptr{Float64}, solution.s)

    x = copy(pointer_to_array(ptr_x, n))
    y = copy(pointer_to_array(ptr_y, m))
    s = copy(pointer_to_array(ptr_s, m))

    SCS_finish(data, p_work)
    return Solution(x, y, s, status)
end
