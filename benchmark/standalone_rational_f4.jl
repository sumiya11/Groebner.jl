
if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @btime gb = Groebner.groebner($system, reduced=false)
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
  60.985 ms (216734 allocations: 21.36 MiB)
noon 6
  265.249 ms (778398 allocations: 85.04 MiB)
noon 7
  2.596 s (4994466 allocations: 494.50 MiB)
eco 5
  1.014 ms (6202 allocations: 2.01 MiB)
eco 7
  11.011 ms (47811 allocations: 6.94 MiB)
katsura 6
  252.392 ms (325027 allocations: 26.47 MiB)
=#
