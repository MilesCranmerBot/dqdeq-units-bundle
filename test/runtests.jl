using Test

using DynamicQuantities
using DiffEqBase
using OrdinaryDiffEqTsit5
using LinearAlgebra

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

        # Keep integration smoke tests focused and lightweight in this bundle:
        # explicit Tsit5 path + unitful callbacks/saveat/reverse-time checks below.
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

    @testset "unitful default_factorize linear solve (stiff-path surrogate)" begin
        # First entry is zero to guard against first-element unit assumptions.
        A = [0.0u"s^-1" 1.0u"s^-1"; -2.0u"s^-1" 3.0u"s^-1"]
        b = [3.0u"m/s", 4.0u"m/s"]

        W = DiffEqBase.default_factorize(A)
        x = W \ b

        @test x isa Vector
        @test eltype(x) <: typeof(1.0u"m")

        # Check Ax ≈ b in unitless space (same physical units on both sides).
        Ax = A * x
        @test ustrip_unit(Ax[1], 1.0u"m/s") ≈ ustrip_unit(b[1], 1.0u"m/s") atol=1e-12 rtol=1e-12
        @test ustrip_unit(Ax[2], 1.0u"m/s") ≈ ustrip_unit(b[2], 1.0u"m/s") atol=1e-12 rtol=1e-12

        # Also exercise the ldiv! path used by Rosenbrock/SDIRK internals.
        y = similar(b)
        LinearAlgebra.ldiv!(y, W, b)
        @test ustrip_unit(y[1], 1.0u"m") ≈ ustrip_unit(x[1], 1.0u"m") atol=1e-12 rtol=1e-12
        @test ustrip_unit(y[2], 1.0u"m") ≈ ustrip_unit(x[2], 1.0u"m") atol=1e-12 rtol=1e-12
    end
end
