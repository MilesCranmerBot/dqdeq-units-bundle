using Test
using DynamicQuantities
using DifferentialEquations

# Helpers
ustrip_unit(x, u) = x / u

@testset "DynamicQuantities × DifferentialEquations: unitful ODEs" begin
    @testset "out-of-place vector state" begin
        u0 = [1.0u"m", 2.0u"m"]
        tspan = (0.0u"s", 1.0u"s")
        f(u, p, t) = u ./ (1u"s")
        prob = ODEProblem(f, u0, tspan)

        function check_final(sol; atol=1e-6, rtol=1e-6)
            @test sol.retcode == ReturnCode.Success
            @test sol.t[end] == 1.0u"s"
            @test ustrip_unit(sol.u[end][1], 1.0u"m") ≈ exp(1) atol=atol rtol=rtol
            @test ustrip_unit(sol.u[end][2], 1.0u"m") ≈ 2 * exp(1) atol=atol rtol=rtol
        end

        # Baseline explicit solver
        sol_tsit = solve(prob, Tsit5())
        check_final(sol_tsit; atol=1e-6, rtol=1e-6)

        # Interpolation should also preserve units
        mid = sol_tsit(0.5u"s")
        @test mid[1] isa typeof(1.0u"m")
        @test ustrip_unit(mid[1], 1.0u"m") ≈ exp(0.5) atol=2e-6 rtol=2e-6

        # A few other explicit/nonstiff solvers
        for alg in (BS3(), Vern7(), DP5())
            sol = solve(prob, alg; reltol=1e-12, abstol=1e-12)
            check_final(sol; atol=2e-9, rtol=2e-9)
        end

        # Stiff solvers (exercise Jacobian/linear-solve paths)
        sol_rodas = solve(prob, Rodas5())
        check_final(sol_rodas; atol=1e-5, rtol=1e-5)

        sol_ros = solve(prob, Rosenbrock23(); reltol=1e-10, abstol=1e-10)
        check_final(sol_ros; atol=2e-7, rtol=2e-7)
    end

    @testset "scalar state" begin
        u0 = 1.0u"m"
        tspan = (0.0u"s", 1.0u"s")
        f(u, p, t) = u / (1u"s")
        prob = ODEProblem(f, u0, tspan)

        sol = solve(prob, Tsit5())
        @test sol.t[end] == 1.0u"s"
        @test ustrip_unit(sol.u[end], 1.0u"m") ≈ exp(1) atol=1e-6 rtol=1e-6
    end

    @testset "reverse-time integration" begin
        u0 = 1.0u"m"
        tspan = (1.0u"s", 0.0u"s")
        f(u, p, t) = u / (1u"s")
        prob = ODEProblem(f, u0, tspan)

        sol = solve(prob, Tsit5())
        @test sol.t[end] == 0.0u"s"
        # u(0) = exp(-1) * u(1)
        @test ustrip_unit(sol.u[end], 1.0u"m") ≈ exp(-1) atol=2e-6 rtol=2e-6
    end

    @testset "saveat is unitful" begin
        u0 = 1.0u"m"
        tspan = (0.0u"s", 1.0u"s")
        f(u, p, t) = u / (1u"s")
        prob = ODEProblem(f, u0, tspan)

        saveat = [0.0u"s", 0.25u"s", 0.5u"s", 0.75u"s", 1.0u"s"]
        sol = solve(prob, Tsit5(); saveat)
        @test sol.t == saveat
        @test ustrip_unit(sol.u[end], 1.0u"m") ≈ exp(1) atol=1e-6 rtol=1e-6
    end

    @testset "callbacks with unitful time" begin
        u0 = 1.0u"m"
        tspan = (0.0u"s", 1.0u"s")
        f(u, p, t) = u / (1u"s")
        prob = ODEProblem(f, u0, tspan)

        hits = Ref(0)
        condition(u, t, integrator) = (t == 0.5u"s")
        affect!(integrator) = (hits[] += 1)
        cb = DiscreteCallback(condition, affect!; save_positions=(false, false))

        sol = solve(
            prob,
            Tsit5();
            callback=cb,
            tstops=[0.5u"s"],
            saveat=[0.0u"s", 0.5u"s", 1.0u"s"],
        )
        @test hits[] == 1
        @test sol.t == [0.0u"s", 0.5u"s", 1.0u"s"]
        @test ustrip_unit(sol.u[end], 1.0u"m") ≈ exp(1) atol=1e-6 rtol=1e-6
    end
end
