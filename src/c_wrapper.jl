export SCS_init, SCS_solve, SCS_finish, SCS_version

macro compat_gc_preserve(args...)
    vars = args[1:end-1]
    body = args[end]
    if VERSION > v"0.7.0-"
        return esc(Expr(:macrocall, Expr(:., :Base, Base.Meta.quot(Symbol("@gc_preserve"))), __source__, args...))
    else
        return esc(body)
    end
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
function SCS_solve(m::Int, n::Int, A::SCSVecOrMatOrSparse, b::Array{Float64},
        c::Array{Float64}, f::Int, l::Int, q::Array{Int}, s::Array{Int},
        ep::Int, ed::Int, p::Array{Float64},
        primal_sol::Vector{Float64}=Float64[],
        dual_sol::Vector{Float64}=Float64[],
        slack::Vector{Float64}=Float64[];
        options...)

    managed_matrix = ManagedSCSMatrix(m, n, A)
    matrix = Ref(SCSMatrix(managed_matrix))
    settings = Ref(SCSSettings(;options...))

    data = Ref(SCSData(m, n, Base.unsafe_convert(Ptr{SCSMatrix}, matrix), pointer(b), pointer(c), Base.unsafe_convert(Ptr{SCSSettings},settings)))

    cone = Ref(SCSCone(f, l, q, s, ep, ed, p))

    info = Ref(SCSInfo())
    p_work = SCS_init(data, cone, info)

    if (:warm_start, true) in options && length(primal_sol) == n && length(dual_sol) == m && length(slack) == m
        x = primal_sol
        y = dual_sol
        s = slack
    else
        x = zeros(n)
        y = zeros(m)
        s = zeros(m)
    end

    solution = SCSSolution(pointer(x), pointer(y), pointer(s))

    @compat_gc_preserve managed_matrix matrix settings b c q s p begin
        status = SCS_solve(p_work, data, cone, solution, info)
    end
    SCS_finish(p_work)
    return Solution(x, y, s, status)

end

# Wrappers for the direct C API.
# Do not call these wrapper methods directly unless you understand the
# use of @gc_preserve in the SCS_solve helper above.

# Take Ref{}s because SCS might modify the structs
function SCS_init(data::Ref{SCSData}, cone::Ref{SCSCone}, info::Ref{SCSInfo})

    p_work = ccall((:scs_init, SCS.scs), Ptr{Void},
        (Ptr{SCSData}, Ptr{SCSCone}, Ptr{SCSInfo}),
        data, cone, info)

    return p_work
end

# solution struct is simple enough, we know it won't be modified by SCS so take by value
function SCS_solve(p_work::Ptr{Void}, data::Ref{SCSData}, cone::Ref{SCSCone}, solution::SCSSolution, info::Ref{SCSInfo})

    status = ccall((:scs_solve, SCS.scs), Int,
        (Ptr{Void}, Ptr{SCSData}, Ptr{SCSCone}, Ref{SCSSolution}, Ptr{SCSInfo}),
        p_work, data, cone, solution, info)

    return status
end

function SCS_finish(p_work::Ptr{Void})
    ccall((:scs_finish, SCS.scs), Void,
        (Ptr{Void}, ),
        p_work)
end

# This one is safe to call
function SCS_version()
    return unsafe_string(ccall((:scs_version, SCS.scs), Cstring, ()))
end
