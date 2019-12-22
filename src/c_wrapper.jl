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
        primal_sol::Vector{Float64}=zeros(n),
        dual_sol::Vector{Float64}=zeros(m),
        slack::Vector{Float64}=zeros(m);
        options...)

    n > 0 || throw(ArgumentError("The number of variables in SCSModel must be greater than 0"))
    m > 0 || throw(ArgumentError("The number of constraints in SCSModel must be greater than 0"))

    ws = (:warm_start=>true) in options
    warmstart_sizes_are_correct = length(primal_sol) == n && length(dual_sol) == length(slack) == m

    if warmstart_sizes_are_correct
        if !ws
            fill!(primal_sol, 0.0)
            fill!(dual_sol, 0.0)
            fill!(slack, 0.0)
        end
    else
        ws && throw(ArgumentError("The provided warmstart doesn't match problem sizes"))
        primal_sol = zeros(n)
        dual_sol = zeros(m)
        slack = zeros(m)
    end

    solution = SCSSolution(pointer(primal_sol), pointer(dual_sol), pointer(slack))

    settings = Base.cconvert(Ref{SCSSettings}, SCSSettings(T; options...))
    managed_matrix = ManagedSCSMatrix(m, n, A)
    data = Base.cconvert(Ref{SCSData}, SCSData(m, n,
        Base.unsafe_convert(Ref{SCSMatrix}, managed_matrix.scsmat), #creates Ptr{SCSMatrix}
        pointer(b), pointer(c),
        Base.unsafe_convert(Ref{SCSSettings}, settings))) #creates Ptr{SCSSettings}
    # unsafe_convert doesn't protect from GC: mmatrix and settings must be GC.@preserved
    cone = Base.cconvert(Ref{SCSCone}, SCSCone(f, l, q, s, ep, ed, p))
    info = Base.cconvert(Ref{SCSInfo}, SCSInfo())

    Base.GC.@preserve managed_matrix settings b c q s p begin
        p_work = SCS_init(T, data, cone, info)
        status = SCS_solve(T, p_work, data, cone, solution, info)
        SCS_finish(T, p_work)
    end

    return Solution(primal_sol, dual_sol, slack, info[], status)
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
                (Ref{SCSData}, Ref{SCSCone}, Ref{SCSInfo}),
                data, cone, info)

            return p_work
        end

        # solution struct is simple enough, we know it won't be modified by SCS so take by value
        function SCS_solve(::Type{$T}, p_work::Ptr{Nothing}, data::Ref{SCSData}, cone::Ref{SCSCone}, solution::SCSSolution, info::Ref{SCSInfo})

            status = ccall((:scs_solve, $lib), Int,
                (Ptr{Nothing}, Ref{SCSData}, Ref{SCSCone}, Ref{SCSSolution}, Ref{SCSInfo}),
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
