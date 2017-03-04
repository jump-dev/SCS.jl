export SCS_init, SCS_solve, SCS_finish, SCS_version


function SCS_version()
    return unsafe_string(ccall((:scs_version, SCS.scs), Cstring, ()))
end

for (T, lib) in zip([SCS.Direct, SCS.Indirect], [SCS.scs, SCS.scsindir])
    @eval begin
        function SCS_init(::Type{$T}, data::SCSData, cone::SCSCone)
            # Initialize the info struct
            info = SCSInfo(0, convert(Int128, 0), convert(Int128, 0), 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

            p_work = ccall((:scs_init, $lib), Ptr{Void},
                (Ptr{SCSData}, Ptr{SCSCone}, Ptr{SCSInfo}),
                &data, &cone, &info)

            return p_work, info
        end

        function SCS_solve(::Type{$T}, p_work::Ptr{Void}, data::SCSData, cone::SCSCone, info::SCSInfo, solution::SCSSolution)
            solution_ptr = pointer([solution])

            info_ptr = pointer([info])
            status = ccall((:scs_solve, $lib), Int,
                (Ptr{Void}, Ptr{SCSData}, Ptr{SCSCone}, Ptr{SCSSolution}, Ptr{SCSInfo}),
                p_work, &data, &cone, solution_ptr, info_ptr)

            solution = unsafe_load(solution_ptr)
            info = unsafe_load(info_ptr)

            return status, solution, info, p_work
        end

        function SCS_solve(::Type{$T}, p_work::Ptr{Void}, data::SCSData, cone::SCSCone, info::SCSInfo)
            solution = SCSSolution(pointer(zeros(data.n)), pointer(zeros(data.m)), pointer(zeros(data.m)))
            return SCS_solve($T, p_work, data, cone, info, solution)
        end

        function SCS_solve(::Type{$T}, data::SCSData, cone::SCSCone, solution::SCSSolution)
            p_work, info = SCS_init($T, data, cone)
            return SCS_solve($T, p_work, data, cone, info, solution)
        end


        function SCS_finish(::Type{$T}, p_work::Ptr{Void})
            ccall((:scs_finish, $lib), Void,
                (Ptr{Void}, ),
                p_work)
        end
    end
end


function SCS_solve(p_work::Ptr{Void}, data::SCSData, cone::SCSCone, info::SCSInfo, solution::SCSSolution)
    solution_ptr = pointer([solution])

    info_ptr = pointer([info])

    status = ccall((:scs_solve, SCS.scs), Int,
        (Ptr{Void}, Ptr{SCSData}, Ptr{SCSCone}, Ptr{SCSSolution}, Ptr{SCSInfo}),
        p_work, &data, &cone, solution_ptr, info_ptr)

    solution = unsafe_load(solution_ptr)
    info = unsafe_load(info_ptr)

    return status, solution, info, p_work
end


function SCS_solve(p_work::Ptr{Void}, data::SCSData, cone::SCSCone, info::SCSInfo)
    solution = SCSSolution(pointer(zeros(data.n)), pointer(zeros(data.m)), pointer(zeros(data.m)))
    return SCS_solve(p_work, data, cone, info, solution)
end


function SCS_solve(data::SCSData, cone::SCSCone, solution::SCSSolution)
    p_work, info = SCS_init(data, cone)
    return SCS_solve(p_work, data, cone, info, solution)
end


function SCS_finish(p_work::Ptr{Void})
    ccall((:scs_finish, SCS.scs), Void,
        (Ptr{Void}, ),
        p_work)
end
