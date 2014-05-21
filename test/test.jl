import SCS

function test_scs()
  SCSInt = SCS.SCSInt
  A = SCS.create_scs_matrix([1.0]')
  temp = [A]
  A_ptr = pointer(temp)
  b = [1.0]
  c = [1.0]
  data = SCS.create_scs_data(m=convert(SCSInt, 1), n=convert(SCSInt, 1), A=A_ptr, b=pointer(b), c=pointer(c))

  q = [convert(SCSInt, 0)]
  s = [convert(SCSInt, 0)]
  one = convert(SCSInt, 1)
  zero = convert(SCSInt, 0)

  cone = SCS.create_scs_cone(one, zero, pointer(q), zero, pointer(s), zero, zero, zero)

  p_work, info = SCS.init(data, cone)
  status, solution, info, p_work = SCS.solve(p_work, data, cone, info)

  @assert status == 1
  SCS.finish(data, p_work)
  return status
end


# TODO: Delete
# function old_test_scs()
#   sol = SCSSolution(C_NULL, C_NULL, C_NULL)

#   zero = convert(SCSInt, 0)
#   info = SCSInfo(zero, convert(Int128, 0), convert(Int128, 0), zero, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#   # info = SCSInfo(zero, C_NULL, zero, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

#   q = [convert(SCSInt, 0)]
#   s = [convert(SCSInt, 0)]
#   one = convert(SCSInt, 1)

#   # TODO: Check what this should be
#   cone = SCSCone(one, zero, pointer(q), zero, pointer(s), zero, zero, zero)

#   sp = sparse([1.0]')

#   values = sp.nzval * 1.0
#   rowval = convert(Array{SCSInt, 1}, sp.rowval - 1)
#   colptr = convert(Array{SCSInt, 1}, sp.colptr - 1)

#   A = SCSMatrix(pointer(values), pointer(rowval), pointer(colptr))

#   temp = [A]
#   p = pointer(temp)
#   c = [1.0]
#   b = [1.0]

#   data = SCSData(one, one, p, pointer(b),  pointer(c), convert(SCSInt, 2500), 1e-3, 1.8, 1e-3, 5.0, 1.5,
#     one, one, zero)

#   info_ptr = [info]
#   info_ptr = pointer(info_ptr)
#   status = ccall((:scs, "../scs/scs.so"), SCSInt, (Ptr{SCSData}, Ptr{SCSCone}, Ptr{SCSSolution}, Ptr{SCSInfo}),
#                   &data, &cone, &sol, info_ptr)
#   return status, sol, info, info_ptr
# end

