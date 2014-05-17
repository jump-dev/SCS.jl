module SCS

export test_scs, SCSMatrix, SCSData, SCSSolution, SCSInfo, SCSCone

immutable SCSMatrix
  values::Ptr{Float64}
  rowval::Ptr{Int32}
  colptr::Ptr{Int32}
end

immutable SCSData
  m::Int32
  n::Int32
  A::Ptr{SCSMatrix}
  b::Ptr{Float64}
  c::Ptr{Float64}
  max_iters::Int32
  eps::Float64
  alpha::Float64
  rho_x::Float64
  scale::Float64
  cg_rate::Float64
  verbose::Int32
  normalize::Int32
  warm_start::Int32
end

immutable SCSSolution
  x::Ptr{None}
  y::Ptr{None}
  s::Ptr{None}
end

immutable SCSInfo
  iter::Int32

  # TODO: Look into this
  status1::Int128
  status2::Int128

  statusVal::Int32
  pobj::Float64
  dobj::Float64
  resPri::Float64
  resDual::Float64
  relGap::Float64
  setupTime::Float64
  solveTime::Float64
end

immutable SCSCone
  f::Int32
  l::Int32
  q::Ptr{Int32}
  qsize::Int32
  s::Ptr{Int32}
  ssize::Int32
  ep::Int32
  ed::Int32
end

function test_scs()
  sol = SCSSolution(C_NULL, C_NULL, C_NULL)

  zero = convert(Int32, 0)
  info = SCSInfo(zero, convert(Int128, 0), convert(Int128, 0), zero, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  # info = SCSInfo(zero, C_NULL, zero, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

  q = [convert(Int32, 0)]
  s = [convert(Int32, 0)]
  one = convert(Int32, 1)

  # TODO: Check what this should be
  cone = SCSCone(one, zero, pointer(q), zero, pointer(s), zero, zero, zero)

  sp = sparse([1.0]')

  values = sp.nzval * 1.0
  rowval = convert(Array{Int32, 1}, sp.rowval - 1)
  colptr = convert(Array{Int32, 1}, sp.colptr - 1)

  A = SCSMatrix(pointer(values), pointer(rowval), pointer(colptr))

  temp = [A]
  p = pointer(temp)
  c = [1.0]
  b = [1.0]

  data = SCSData(one, one, p, pointer(b),  pointer(c), convert(Int32, 2500), 1e-3, 1.8, 1e-3, 5.0, 1.5,
    one, one, zero)

  status = ccall((:scs, "../scs/scs.so"), Int32, (Ptr{SCSData}, Ptr{SCSCone}, Ptr{SCSSolution}, Ptr{SCSInfo}),
                  &data, &cone, &sol, &info)

end


end # module
