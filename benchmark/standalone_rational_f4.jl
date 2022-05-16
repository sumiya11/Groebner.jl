
if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BenchmarkTools.DEFAULT_PARAMETERS.samples = 3

function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @btime gb = Groebner.groebner($system, reduced=false)
    # gb = Groebner.groebner(system, reduced=false)
end

function benchmark_system_my_certify(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @btime gb = Groebner.groebner($system, reduced=false, linalg=:prob)
    # gb = Groebner.groebner(system, reduced=false)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 11", Groebner.rootn(11, ground=ground)),
        ("noon 7"    ,Groebner.noonn(7, ground=ground)),
        ("eco 10"    ,Groebner.eco10(ground=ground)),
        ("katsura 8"    ,Groebner.katsuran(8, ground=ground))
        ]

    println()
    for (name, system) in systems
        println("$name")
        # benchmark_system_my(system)
        benchmark_system_my_certify(system)
    end
end

ground = AbstractAlgebra.QQ
run_f4_ff_degrevlex_benchmarks(ground)

#=
cyclic 11
  50.711 ms (191971 allocations: 23.54 MiB)
noon 7
  910.447 ms (2087052 allocations: 298.79 MiB)
eco 10
  1.484 s (1994909 allocations: 417.28 MiB)
katsura 8
  4.509 s (2982135 allocations: 499.90 MiB)
=#
