using SCS
using Base.Test

# Solve a trivial problem
A = reshape([1.0],(1,1))
solution = SCS_solve(1, 1, A, [1.0], [1.0], 1, 0, Int[], Int[], 0, 0, Float64[]);
@assert solution.ret_val == 1




