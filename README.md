# DynamicQuantities × DifferentialEquations: unitful ODE repro bundle

This repo vendors snapshots of a small set of SciML/DynamicQuantities packages.

Commits:
- **base**: unpatched upstream snapshots
- **patch**: minimal patch stack to make DynamicQuantities unitful time/state work through `DifferentialEquations.solve`, including a stiff Rosenbrock method.

## Run (Julia 1.12)

```bash
julia --version
julia --project=. -e "import Pkg; Pkg.instantiate()"
julia --project=. repro_vector_stiff.jl
```

Expected output includes both Tsit5 and Rodas5 results with unitful `tend = 1.0 s`.
