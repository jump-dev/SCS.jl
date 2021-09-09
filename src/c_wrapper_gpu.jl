if haskey(ENV, "JULIA_SCS_LIBRARY_PATH")
    @isdefined(libscsgpuindir) && push!(available_solvers, GpuIndirectSolver)
else
    import SCS_GPU_jll
    const gpuindirect = SCS_GPU_jll.libscsgpuindir
    push!(available_solvers, GpuIndirectSolver)
end
for linear_solver in (GpuIndirectSolver,)
    # lib = gpuindirect # clib(linear_solver)
    T = scsint_t(linear_solver)
    @eval begin
        function SCS_set_default_settings(
            ::Type{$linear_solver},
            data::SCSData{$T},
        )
            return ccall(
                (:scs_set_default_settings, $(clib(linear_solver))),
                Cvoid,
                (Ref{SCSData{$T}},),
                data,
            )
        end

        # data and cone are const in :scs_init
        function SCS_init(
            ::Type{$linear_solver},
            data::SCSData{$T},
            cone::SCSCone{$T},
            info_ref::Ref{SCSInfo{$T}},
        )
            p_work = ccall(
                (:scs_init, $(clib(linear_solver))),
                Ptr{Cvoid},
                (Ref{SCSData{$T}}, Ref{SCSCone{$T}}, Ref{SCSInfo{$T}}),
                data,
                cone,
                info_ref,
            )

            return p_work
        end

        # data and cone are const in :scs_solve
        # solution struct contains only `Ptr`s so passing by value
        function SCS_solve(
            ::Type{$linear_solver},
            p_work::Ptr{Nothing},
            data::SCSData{$T},
            cone::SCSCone{$T},
            solution::SCSSolution,
            info_ref::Ref{SCSInfo{$T}},
        )
            status = ccall(
                (:scs_solve, $(clib(linear_solver))),
                $T,
                (
                    Ptr{Cvoid},
                    Ref{SCSData{$T}},
                    Ref{SCSCone{$T}},
                    Ref{SCSSolution},
                    Ref{SCSInfo{$T}},
                ),
                p_work,
                data,
                cone,
                solution,
                info_ref,
            )

            return status
        end

        function SCS_finish(::Type{$linear_solver}, p_work::Ptr{Nothing})
            return ccall(
                (:scs_finish, $(clib(linear_solver))),
                Cvoid,
                (Ptr{Cvoid},),
                p_work,
            )
        end
    end
end
