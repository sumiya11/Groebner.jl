using Pkg
Pkg.status()
commit = ARGS[1]
if commit == "master"
    Pkg.add(url="https://github.com/sumiya11/Groebner.jl")
else
    Pkg.add(url="https://github.com/sumiya11/Groebner.jl", rev=commit)
end
# Pkg.update("Groebner")

include("../run_benchmarks.jl")

dump_results((@__DIR__) * "/results", "stable")
