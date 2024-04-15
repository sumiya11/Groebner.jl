import Pkg
Pkg.activate(@__DIR__)
Pkg.develop(path=(@__DIR__) * "/../../../")
Pkg.instantiate()
Pkg.status()

id = ARGS[1]

include("../run_benchmarks.jl")

dump_results((@__DIR__) * "/results_$id", "nightly")
