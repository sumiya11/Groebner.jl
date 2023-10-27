import Pkg
Pkg.develop(path=(@__DIR__) * "/../../../")

include("../run_benchmarks.jl")

dump_results((@__DIR__) * "/results", "nightly")
