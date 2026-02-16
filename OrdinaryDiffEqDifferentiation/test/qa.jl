using Pkg
Pkg.add("Aqua")

using OrdinaryDiffEqDifferentiation
using Aqua

@testset "Aqua" begin
    # Aqua's `persistent_tasks` check runs a separate Julia process that precompiles a
    # wrapper package. On Julia 1.12 we've observed occasional crashes in that
    # subprocess (e.g. `done.log was not created, but precompilation exited`), which
    # makes CI flaky/unstable. Disable it on Julia ≥ 1.12 until the upstream issue is
    # resolved.
    persistent_tasks = VERSION < v"1.12.0-"

    Aqua.test_all(
        OrdinaryDiffEqDifferentiation;
        piracies = false,
        ambiguities = false,
        persistent_tasks,
    )
end
