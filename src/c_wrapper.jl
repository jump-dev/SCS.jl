export SCS_init, SCS_solve, SCS_finish, SCS_version

# SCS solves a problem of the form
# minimize        c' * x
# subject to      A * x + s = b
#                 s in K
# where K is a product cone of
# zero cones,
# linear cones { x | x >= 0 },
# second-order cones (SOC) { (t,x) | ||x||_2 <= t },
# semi-definite cones (SDC) { X | X psd }, and
# exponential cones {(x,y,z) | y e^(x/y) <= z, y>0 }.
#
#
# Description of input argments:
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
# Returns a Solution object.
function SCS_solve(T::Union{Type{Direct}, Type{Indirect}},
        m::Int, n::Int, A::SCSVecOrMatOrSparse, b::Array{Float64},
        c::Array{Float64}, f::Int, l::Int, q::Array{Int}, s::Array{Int},
        ep::Int, ed::Int, p::Array{Float64},
        primal_sol::Vector{Float64}=Float64[],
        dual_sol::Vector{Float64}=Float64[],
        slack::Vector{Float64}=Float64[];
        options...)

    n > 0 || throw(ArgumentError("The number of variables in SCSModel must be greater than 0"))
    m > 0 || throw(ArgumentError("The number of constraints in SCSModel must be greater than 0"))

    managed_matrix = ManagedSCSMatrix(m, n, A)
    matrix = Ref(SCSMatrix(managed_matrix))
    settings = Ref(SCSSettings(T; options...))
    data = Ref(SCSData(m, n, Base.unsafe_convert(Ptr{SCSMatrix}, matrix), pointer(b), pointer(c), Base.unsafe_convert(Ptr{SCSSettings},settings)))

    cone = Ref(SCSCone(f, l, q, s, ep, ed, p))
    info = Ref(SCSInfo())

    ws = (:warm_start=>true) in options

    if ws && length(primal_sol) == n && length(dual_sol) == m && length(slack) == m
        x = primal_sol
        y = dual_sol
        s = slack
    else
        x = zeros(n)
        y = zeros(m)
        s = zeros(m)
    end
    solution = SCSSolution(pointer(x), pointer(y), pointer(s))

    Base.GC.@preserve managed_matrix matrix settings b c q s p begin
        p_work = SCS_init(T, data, cone, info)
        status = SCS_solve(T, p_work, data, cone, solution, info)
        SCS_finish(T, p_work)
    end

    return Solution(x, y, s, info[], status)

end

# Wrappers for the direct C API.
# Do not call these wrapper methods directly unless you understand the
# use of @gc_preserve in the SCS_solve helper above.

# Take Ref{}s because SCS might modify the structs
for (T, lib) in zip([SCS.Direct, SCS.Indirect], [SCS.direct, SCS.indirect])
    @eval begin

        function SCS_set_default_settings(::Type{$T}, data::Ref{SCSData})
            ccall((:scs_set_default_settings, $lib), Nothing, (Ref{SCSData}, ), data)
        end

        function SCS_init(::Type{$T}, data::Ref{SCSData}, cone::Ref{SCSCone}, info::Ref{SCSInfo})

            p_work = ccall((:scs_init, $lib), Ptr{Nothing},
                (Ptr{SCSData}, Ptr{SCSCone}, Ptr{SCSInfo}),
                data, cone, info)

            return p_work
        end

        # solution struct is simple enough, we know it won't be modified by SCS so take by value
        function SCS_solve(::Type{$T}, p_work::Ptr{Nothing}, data::Ref{SCSData}, cone::Ref{SCSCone}, solution::SCSSolution, info::Ref{SCSInfo})

            status = ccall((:scs_solve, $lib), Int,
                (Ptr{Nothing}, Ptr{SCSData}, Ptr{SCSCone}, Ref{SCSSolution}, Ptr{SCSInfo}),
                p_work, data, cone, solution, info)

            return status
        end

        function SCS_finish(::Type{$T}, p_work::Ptr{Nothing})
            ccall((:scs_finish, $lib), Nothing,
                (Ptr{Nothing}, ),
                p_work)
        end
    end
end

# This one is safe to call
function SCS_version()
    return unsafe_string(ccall((:scs_version, SCS.direct), Cstring, ()))
end
