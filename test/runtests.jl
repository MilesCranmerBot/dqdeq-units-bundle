using Test
using DynamicQuantities
using DifferentialEquations

@testset "DynamicQuantities × DifferentialEquations (unitful ODE)" begin
    u0 = [1.0u"m", 2.0u"m"]
    tspan = (0.0u"s", 1.0u"s")
    f(u,p,t) = u ./ (1u"s")

    prob = ODEProblem(f, u0, tspan)

    sol_tsit = solve(prob, Tsit5())
    @test sol_tsit.t[end] == 1.0u"s"
    @test sol_tsit.u[end][1] / (1.0u"m") ≈ exp(1) atol=1e-6 rtol=1e-6
    @test sol_tsit.u[end][2] / (1.0u"m") ≈ 2 * exp(1) atol=1e-6 rtol=1e-6

    sol_rodas = solve(prob, Rodas5())
    @test sol_rodas.t[end] == 1.0u"s"
    @test sol_rodas.u[end][1] / (1.0u"m") ≈ exp(1) atol=1e-5 rtol=1e-5
    @test sol_rodas.u[end][2] / (1.0u"m") ≈ 2 * exp(1) atol=1e-5 rtol=1e-5
end
