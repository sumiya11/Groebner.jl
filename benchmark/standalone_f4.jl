
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
68.460 ms (184276 allocations: 29.25 MiB)
cyclic 13
255.064 ms (459892 allocations: 74.68 MiB)
katsura 9
344.615 ms (82539 allocations: 22.35 MiB)
katsura 10
2.725 s (250078 allocations: 74.80 MiB)
noon 6
36.954 ms (101244 allocations: 21.19 MiB)
noon 7
289.916 ms (494536 allocations: 100.91 MiB)
=#
