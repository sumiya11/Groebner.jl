
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

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 7", Groebner.cyclicn(7, ground=ground)),
        ("cyclic 8", Groebner.cyclicn(8, ground=ground)),
        ("katsura 10",Groebner.katsura10(ground=ground)),
        ("katsura 11",Groebner.katsura11(ground=ground)),
        ("eco 12",Groebner.eco12(ground=ground)),
        ("eco 13",Groebner.eco13(ground=ground)),
        ("noon 7"    ,Groebner.noonn(7, ground=ground)),
        ("noon 8"    ,Groebner.noonn(8, ground=ground))
    ]

    println("Running Groebner.jl benchmarks")
    for (name, system) in systems
        println("==================")
        println("System $name:")
        benchmark_system_groebner(system)
    end
end

ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)
