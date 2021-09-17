using LinearAlgebra

args = let
    #= the same as
    using JuMP, LinearAlgebra
    m = let
        m = Model()
        @variable m x[1:5] >=0;
        @constraints m begin
            x[1] + x[2] <= 5
            x[2] + x[5] <= 3
            x[3] + x[4] + x[5] <= 9
        end
        c = - [3.0, 4.0, 4.0, 9.0, 5.0]
        @objective m Min dot(c, x)
        m
    end
    =#

    m = 8
    n = 5
    A = [
        1.0 1.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0 1.0
        0.0 0.0 1.0 1.0 1.0
        -1.0 0.0 0.0 0.0 0.0
        0.0 -1.0 0.0 0.0 0.0
        0.0 0.0 -1.0 0.0 0.0
        0.0 0.0 0.0 -1.0 0.0
        0.0 0.0 0.0 0.0 -1.0
    ]
    P = zeros(n, n)
    b = [5.0, 3.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    c = -[3.0, 4.0, 4.0, 9.0, 5.0]
    z = 0
    l = 8
    bu = Float64[]
    bl = Float64[]
    q = Int64[]
    s = Int64[]
    ep = 0
    ed = 0
    p = Float64[]

    (
        m = m,
        n = n,
        A = A,
        P = P,
        b = b,
        c = c,
        z = z,
        l = l,
        bu = bu,
        bl = bl,
        q = q,
        s = s,
        ep = ep,
        ed = ed,
        p = p,
    )
end

solution = SCS.scs_solve(SCS.DirectSolver, args...)

@test isapprox(dot(solution.x, args.c), -99.0, rtol = 1e-6)
@test !isapprox(dot(solution.x, args.c), -99.0, rtol = 1e-7)

solution =
    SCS.scs_solve(SCS.DirectSolver, args..., eps_abs = 1e-12, eps_rel = 1e-10)
@test isapprox(dot(solution.x, args.c), -99.0, rtol = 1e-12)

solution = SCS.scs_solve(
    SCS.DirectSolver,
    args...,
    solution.x,
    solution.y,
    solution.s,
    max_iters = 2,
    warm_start = true,
    eps_abs = 1e-12,
    eps_rel = 1e-12,
)
@test isapprox(dot(solution.x, args.c), -99.0, rtol = 1e-12)

@test_throws ArgumentError SCS.scs_solve(SCS.DirectSolver, args..., eps = 1e-12)

err = try
    SCS.scs_solve(SCS.DirectSolver, args..., eps_abs = 1e-12, eps = 1e-12)
catch ex
    ex
end
@test err.msg == "Unrecognized option passed to SCS solver: eps"
