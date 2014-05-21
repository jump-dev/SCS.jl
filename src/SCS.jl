module SCS

export SCSInt, create_scs_matrix, create_scs_data, create_scs_cone, init, solve, finish

# TODO: Actually fix this and make it work for windows
@windows_only SCSInt = Int64
@unix_only SCSInt = Int32


immutable SCSMatrix
  values::Ptr{Cdouble}
  rowval::Ptr{SCSInt}
  colptr::Ptr{SCSInt}
end


immutable SCSData
  m::SCSInt
  n::SCSInt
  A::Ptr{SCSMatrix}
  b::Ptr{Cdouble}
  c::Ptr{Cdouble}
  max_iters::SCSInt
  eps::Cdouble
  alpha::Cdouble
  rho_x::Cdouble
  scale::Cdouble
  cg_rate::Cdouble
  verbose::SCSInt
  normalize::SCSInt
  warm_start::SCSInt
end


immutable SCSSolution
  x::Ptr{None}
  y::Ptr{None}
  s::Ptr{None}
end


immutable SCSInfo
  iter::SCSInt
  # We need to allocate 32 bytes for a character string, so we allocate 256 bits
  # of integer instead
  # TODO: Find a better way to do this
  status1::Int128
  status2::Int128

  statusVal::SCSInt
  pobj::Cdouble
  dobj::Cdouble
  resPri::Cdouble
  resDual::Cdouble
  relGap::Cdouble
  setupTime::Cdouble
  solveTime::Cdouble
end


immutable SCSCone
  f::SCSInt
  l::SCSInt
  q::Ptr{SCSInt}
  qsize::SCSInt
  s::Ptr{SCSInt}
  ssize::SCSInt
  ep::SCSInt
  ed::SCSInt
end


immutable SCSWork
  u::Ptr{Cdouble}
  v::Ptr{Cdouble}
  u_t::Ptr{Cdouble}
  u_prev::Ptr{Cdouble}
  h::Ptr{Cdouble}
  g::Ptr{Cdouble}
  pr::Ptr{Cdouble}
  dr::Ptr{Cdouble}
  gTh::Cdouble
  sc_b::Cdouble
  sc_c::Cdouble
  nm_b::Cdouble
  nm_c::Cdouble
  meanNormRowA::Cdouble
  D::Ptr{Cdouble}
  E::Ptr{Cdouble}
  p::Ptr{Void}
end


function init(data::SCSData, cone::SCSCone)
  # Initialize the info struct
  zero = convert(SCSInt, 0)
  info = SCSInfo(zero, convert(Int128, 0), convert(Int128, 0), zero, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

  p_work = ccall((:scs_init, "../scs/scs.so"), Ptr{SCSWork},
      (Ptr{SCSData}, Ptr{SCSCone}, Ptr{SCSInfo}),
      &data, &cone, &info)

  return p_work, info
end


function solve(p_work::Ptr{SCSWork}, data::SCSData, cone::SCSCone, info::SCSInfo)
  # Initialize the solution struct. Note that the pointers can be null since SCS will take care
  # of it for us
  solution = SCSSolution(C_NULL, C_NULL, C_NULL)
  solution_ptr = pointer([solution])

  info_ptr = pointer([info])

  status = ccall((:scs_solve, "../scs/scs.so"), SCSInt,
      (Ptr{SCSWork}, Ptr{SCSData}, Ptr{SCSCone}, Ptr{SCSSolution}, Ptr{SCSInfo}),
      p_work, &data, &cone, solution_ptr, info_ptr)

  solution = unsafe_load(solution_ptr)
  info = unsafe_load(info_ptr)


  return status, solution, info, p_work
end


function finish(data::SCSData, p_work::Ptr{SCSWork})
  ccall((:scs_finish, "../scs/scs.so"), Void,
      (Ptr{SCSData}, Ptr{SCSWork}),
      &data, p_work)
end


# All of these will be keyword arguments
function create_scs_data(;m::SCSInt=nothing, n::SCSInt=nothing, A::Ptr{SCSMatrix}=nothing,
    b::Ptr{Cdouble}=nothing,  c::Ptr{Cdouble}=nothing, max_iters=convert(SCSInt, 2500)::SCSInt,
    eps=convert(Cdouble, 1e-3)::Cdouble, alpha=convert(Cdouble, 1.8)::Cdouble,
    rho_x=convert(Cdouble, 1e-3)::Cdouble, scale=convert(Cdouble, 5.0)::Cdouble,
    cg_rate=convert(Cdouble, 1.5)::Cdouble, verbose=convert(SCSInt, 1)::SCSInt,
    normalize=convert(SCSInt, 1)::SCSInt, warm_start=convert(SCSInt, 0)::SCSInt)

  data = SCSData(m, n, A, b, c, max_iters, eps, alpha, rho_x, scale, cg_rate, verbose, normalize, warm_start)
  return data
end


function create_scs_data(m::SCSInt, n::SCSInt, A::Ptr{SCSMatrix}, b::Ptr{Cdouble}, c::Ptr{Cdouble})
  return create_scs_data(m=m, n=n, A=A, b=b, c=c)
end


function create_scs_matrix(A)
  A_sparse = sparse(A)

  values = A_sparse.nzval * 1.0
  rowval = convert(Array{SCSInt, 1}, A_sparse.rowval - 1)
  colptr = convert(Array{SCSInt, 1}, A_sparse.colptr - 1)

  return SCSMatrix(pointer(values), pointer(rowval), pointer(colptr))
end


function create_scs_cone(f::SCSInt, l::SCSInt, q::Ptr{SCSInt}, qsize::SCSInt, s::Ptr{SCSInt},
  ssize::SCSInt, ep::SCSInt, ed::SCSInt)
  return SCSCone(f, l, q, qsize, s, ssize, ep, ed)
end


end # module
