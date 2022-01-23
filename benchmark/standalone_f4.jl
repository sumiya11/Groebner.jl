
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
    println("length = $(length(gb))")
    println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @btime gb = Groebner.groebner($system, reduced=false)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 12", Groebner.rootn(12, ground=ground)),
        ("cyclic 13", Groebner.rootn(13, ground=ground)),
        ("katsura 9", Groebner.katsura9(ground=ground)),
        ("katsura 10",Groebner.katsura10(ground=ground)),
        ("noon 6"    ,Groebner.noonn(6, ground=ground)),
        ("noon 7"    ,Groebner.noonn(7, ground=ground))
    ]

    for (name, system) in systems
        println("$name")
        benchmark_system_my(system)
    end
end

ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)

#=
cyclic 12
67.473 ms (184149 allocations: 28.37 MiB)
cyclic 13
222.092 ms (459739 allocations: 70.34 MiB)
katsura 9
320.387 ms (82344 allocations: 21.07 MiB)
katsura 10
2.562 s (249743 allocations: 71.48 MiB)
noon 6
36.682 ms (100998 allocations: 21.99 MiB)
noon 7
315.354 ms (493959 allocations: 99.07 MiB)
=#
