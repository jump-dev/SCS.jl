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

    data = SCSData(m, n, Base.unsafe_convert(Ptr{SCSMatrix}, matrix), pointer(b), pointer(c), Base.unsafe_convert(Ptr{SCSSettings}, settings))

    cone = SCSCone(f, l, q, s, ep, ed, p)

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

    status = SCS_solve(p_work, data, cone, solution, info)
    SCS_finish(p_work)
    return Solution(x, y, s, status)

end
