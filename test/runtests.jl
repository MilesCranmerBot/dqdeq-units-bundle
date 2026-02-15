using Test
using DynamicQuantities
using DifferentialEquations

# Helpers
ustrip_unit(x, u) = x / u

@testset "DynamicQuantities × DifferentialEquations: unitful ODEs" begin
    @testset "out-of-place vector state (Tsit5 + Rodas5)" begin
        u0 = [1.0u"m", 2.0u"m"]
        tspan = (0.0u"s", 1.0u"s")
        f(u, p, t) = u ./ (1u"s")
        prob = ODEProblem(f, u0, tspan)

        sol_tsit = solve(prob, Tsit5())
        @test sol_tsit.t[end] == 1.0u"s"
        @test ustrip_unit(sol_tsit.u[end][1], 1.0u"m") ≈ exp(1) atol=1e-6 rtol=1e-6
        @test ustrip_unit(sol_tsit.u[end][2], 1.0u"m") ≈ 2 * exp(1) atol=1e-6 rtol=1e-6

        # Interpolation should also preserve units
        mid = sol_tsit(0.5u"s")
        @test mid[1] isa typeof(1.0u"m")
        @test ustrip_unit(mid[1], 1.0u"m") ≈ exp(0.5) atol=2e-6 rtol=2e-6

        sol_rodas = solve(prob, Rodas5())
        @test sol_rodas.t[end] == 1.0u"s"
        @test ustrip_unit(sol_rodas.u[end][1], 1.0u"m") ≈ exp(1) atol=1e-5 rtol=1e-5
        @test ustrip_unit(sol_rodas.u[end][2], 1.0u"m") ≈ 2 * exp(1) atol=1e-5 rtol=1e-5
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
