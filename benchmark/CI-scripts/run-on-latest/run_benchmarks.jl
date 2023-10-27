import Pkg
Pkg.add("Groebner")

include("../run_benchmarks.jl")

dump_results((@__DIR__) * "/results", "latest")
