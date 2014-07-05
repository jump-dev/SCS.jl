export SCSMatrix, SCSData, SCSSolution, SCSInfo, SCSCone, SCSWork, SCSVecOrMatOrSparse


SCSVecOrMatOrSparse = Union(VecOrMat, SparseMatrixCSC{Float64,Int64})


immutable SCSMatrix
  values::Ptr{Cdouble}
  rowval::Ptr{Clong}
  colptr::Ptr{Clong}
end


immutable SCSData
  # A has m rows, n cols
  m::Clong
  n::Clong
  A::Ptr{SCSMatrix}
  # b is of size m, c is of size n
  b::Ptr{Cdouble}
  c::Ptr{Cdouble}
  # max_iters to take
  max_iters::Clong
  # convergence tolerance
  eps::Cdouble
  # relaxation parameter
  alpha::Cdouble
  # x equality constraint scaling
  rho_x::Cdouble
  # if normalized, rescales by this factor
  scale::Cdouble
  # for indirect, tolerance goes down like (1/iter)^CG_RATE: 1.5
  cg_rate::Cdouble
  verbose::Clong
  # 0 or 1
  normalize::Clong
  # 0 or 1
  warm_start::Clong
end


immutable SCSSolution
  x::Ptr{None}
  y::Ptr{None}
  s::Ptr{None}
end


immutable SCSInfo
  iter::Clong
  # We need to allocate 32 bytes for a character string, so we allocate 256 bits
  # of integer instead
  # TODO: Find a better way to do this
  status1::Int128
  status2::Int128

  statusVal::Clong
  pobj::Cdouble
  dobj::Cdouble
  resPri::Cdouble
  resDual::Cdouble
  relGap::Cdouble
  setupTime::Cdouble
  solveTime::Cdouble
end


immutable SCSCone
  f::Clong # number of linear equality constraints
  l::Clong # length of LP cone
  q::Ptr{Clong} # array of second-order cone constraints
  qsize::Clong # length of SOC array
  s::Ptr{Clong} # array of semi-definite constraints
  ssize::Clong # length of SD array
  ep::Clong # number of primal exponential cone triples
  ed::Clong # number of dual exponential cone triples
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


type Solution
  x::Array{Float64, 1}
  y::Array{Float64, 1}
  s::Array{Float64, 1}
  status::ASCIIString
  ret_val::Int64

  const status_map = {
    1 => "solved",
    -2 => "primal infeasible, dual unbounded",
    -1 => "primal unbounded, dual infeasible",
    -3 => "indeterminate",
    -4 => "failure"
  }

  function Solution(x::Array{Float64, 1}, y::Array{Float64, 1}, s::Array{Float64, 1}, ret_val::Int64)
    if haskey(status_map, ret_val)
      return new(x, y, s, status_map[ret_val], ret_val)
    else
      return new(x, y, s, "unknown problem in solver", ret_val)
    end
  end
end
