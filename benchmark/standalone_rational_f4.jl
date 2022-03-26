
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

    @btime gb = Groebner.groebner($system, linalg=:prob, reduced=false)
    # gb = Groebner.groebner(system, reduced=false)
end

function benchmark_system_my_certify(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @btime gb = Groebner.groebner($system, certify=true, reduced=false)
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
        # benchmark_system_my(system)
        benchmark_system_my_certify(system)
    end
end

ground = AbstractAlgebra.QQ
run_f4_ff_degrevlex_benchmarks(ground)

#=
cyclic 11
  51.766 ms (200404 allocations: 20.83 MiB)
noon 6
  159.142 ms (413056 allocations: 61.13 MiB)
noon 7
  1.189 s (2106605 allocations: 285.85 MiB)
eco 5
  563.600 μs (3556 allocations: 1.80 MiB)
eco 7
  6.228 ms (26509 allocations: 5.09 MiB)
katsura 6
  99.311 ms (94019 allocations: 13.37 MiB)
=#
