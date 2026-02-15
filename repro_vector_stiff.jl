import Pkg
Pkg.activate(@__DIR__)

using DynamicQuantities, DifferentialEquations

u0 = [1.0u"m", 2.0u"m"]
tspan = (0.0u"s", 1.0u"s")
f(u,p,t) = u ./ (1u"s")
prob = ODEProblem(f, u0, tspan)

sol1 = solve(prob, Tsit5())
println("tsit5_end=" * string(sol1.u[end]))
println("tsit5_tend=" * string(sol1.t[end]))

sol2 = solve(prob, Rodas5())
println("rodas5_end=" * string(sol2.u[end]))
println("rodas5_tend=" * string(sol2.t[end]))
