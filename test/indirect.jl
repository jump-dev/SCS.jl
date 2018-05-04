# Solve a trivial problem
A = reshape([1.0],(1,1))
solution = SCS_solve(SCS.Indirect, 1, 1, A, [1.0], [1.0], 1, 0, Int[], Int[], 0, 0, Float64[]);
@test solution.ret_val == 1

feasible_basic_conic(SCS.Indirect)

feasible_exponential_conic(SCS.Indirect);

feasible_sdp_conic(SCS.Indirect)

feasible_pow_conic(SCS.Indirect)
