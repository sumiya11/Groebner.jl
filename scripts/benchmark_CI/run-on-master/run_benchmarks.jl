using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
Pkg.status()

Pkg.add(url="https://github.com/sumiya11/Groebner.jl")

id = ARGS[1]

include("../common.jl")
include("../run_benchmarks.jl")

dump_results((@__DIR__) * "/results_$id", "master")
