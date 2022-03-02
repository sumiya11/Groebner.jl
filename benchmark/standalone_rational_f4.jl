
if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10

function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @btime gb = Groebner.groebner($system, reduced=false)
    # gb = Groebner.groebner(system, reduced=false)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 11", Groebner.rootn(11, ground=ground)),
        ("noon 6"    ,Groebner.noonn(6, ground=ground)),
        ("noon 7"    ,Groebner.noonn(7, ground=ground)),
        ("eco 5"    ,Groebner.eco5(ground=ground)),
        ("eco 7"    ,Groebner.eco7(ground=ground)),
        ("katsura 6"    ,Groebner.katsura6(ground=ground)),
    ]

    println()
    for (name, system) in systems
        println("$name")
        benchmark_system_my(system)
    end
end

ground = AbstractAlgebra.QQ
run_f4_ff_degrevlex_benchmarks(ground)

#=
cyclic 11
  50.663 ms (200419 allocations: 21.05 MiB)
noon 6
  198.116 ms (574503 allocations: 81.15 MiB)
noon 7
  1.847 s (3571431 allocations: 467.36 MiB)
eco 5
  878.300 μs (5467 allocations: 2.00 MiB)
eco 7
  9.300 ms (38332 allocations: 6.76 MiB)
katsura 6
  194.193 ms (211652 allocations: 24.31 MiB)
=#
