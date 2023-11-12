import Pkg
Pkg.add("Groebner")
Pkg.update("Groebner")

include("../run_benchmarks.jl")

dump_results((@__DIR__) * "/results", "stable")
