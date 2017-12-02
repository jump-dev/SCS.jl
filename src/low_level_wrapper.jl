export SCS_init, SCS_solve, SCS_finish, SCS_version


function SCS_version()
    return unsafe_string(ccall((:scs_version, SCS.scs), Cstring, ()))
end


function SCS_init(data::SCSData, cone::SCSCone, info::Ref{SCSInfo})

    p_work = ccall((:scs_init, SCS.scs), Ptr{Void},
        (Ref{SCSData}, Ref{SCSCone}, Ref{SCSInfo}),
        data, cone, info)

    return p_work
end


function SCS_solve(p_work::Ptr{Void}, data::SCSData, cone::SCSCone, solution::SCSSolution, info::Ref{SCSInfo})

    status = ccall((:scs_solve, SCS.scs), Int,
        (Ptr{Void}, Ref{SCSData}, Ref{SCSCone}, Ref{SCSSolution}, Ref{SCSInfo}),
        p_work, data, cone, solution, info)

    return status
end


function SCS_finish(p_work::Ptr{Void})
    ccall((:scs_finish, SCS.scs), Void,
        (Ptr{Void}, ),
        p_work)
end
