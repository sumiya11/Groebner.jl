
# Load all required packages
import Pkg
Pkg.add("Groebner")
Pkg.add("AbstractAlgebra")
Pkg.add("BenchmarkTools")

# Import the packages
using Groebner
import AbstractAlgebra
using BenchmarkTools
using Logging
include((@__DIR__)*"/generate/benchmark_systems.jl")
global_logger(ConsoleLogger(stderr, Logging.Error))

# Benchmark the given system
function benchmark_system_groebner(system)
    system = Groebner.change_ordering(system, :degrevlex)
    # compile:
    Groebner.groebner([system[1]])
    gb = Groebner.groebner(system)
    # run the actual benchmark:
    @btime gb = Groebner.groebner($system)
end

function run_f4_ff_degrevlex_benchmarks(flag)
    if flag
        ground = AbstractAlgebra.GF(2^31-1)
        systems = benchmark_systems_ff(ground)
    else
        ground = AbstractAlgebra.QQ
        systems = benchmark_systems_qq(ground)
    end

    println("Running Groebner.jl benchmarks over $ground")
    for (name, system) in systems
        println("==================")
        println("System $name:")
        benchmark_system_groebner(system)
    end
end

run_f4_ff_degrevlex_benchmarks(true)
run_f4_ff_degrevlex_benchmarks(false)
