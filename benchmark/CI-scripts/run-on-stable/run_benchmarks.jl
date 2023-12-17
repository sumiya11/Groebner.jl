import Pkg
Pkg.status()
Pkg.add(url="https://github.com/sumiya11/Groebner.jl")
# Pkg.update("Groebner")

include("../run_benchmarks.jl")

dump_results((@__DIR__) * "/results", "stable")
