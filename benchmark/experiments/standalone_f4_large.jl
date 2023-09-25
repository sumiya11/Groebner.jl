
using Groebner

import AbstractAlgebra
using BenchmarkTools
using Logging
# global_logger(ConsoleLogger(stderr, Logging.Error))

Groebner.logging_enabled() = false
Groebner.invariants_enabled() = false

function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)

    Groebner.groebner(system)

    t = Inf
    for i in 1:3
        timing = @timed Groebner.groebner(system)
        t = min(t, timing.time)
    end
    t
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 7", Groebner.cyclicn(7, ground=ground)),
        ("cyclic 8", Groebner.cyclicn(8, ground=ground)),
        ("cyclic 9", Groebner.cyclicn(9, ground=ground)),
        ("katsura 9", Groebner.katsuran(9, ground=ground)),
        ("katsura 10", Groebner.katsuran(10, ground=ground)),
        ("katsura 11", Groebner.katsuran(11, ground=ground)),
        ("katsura 12", Groebner.katsuran(12, ground=ground)),
        ("noon 7", Groebner.noonn(7, ground=ground)),
        ("noon 8", Groebner.noonn(8, ground=ground)),
        ("noon 9", Groebner.noonn(9, ground=ground)),
        ("eco 11", Groebner.eco11(ground=ground)),
        ("eco 12", Groebner.eco12(ground=ground)),
        ("eco 12", Groebner.eco13(ground=ground))
    ]

    for (name, system) in systems
        print("$name -- ")
        t = benchmark_system_my(system)
        println("$t s")
    end
end

ground = AbstractAlgebra.GF(2^30 + 3)
run_f4_ff_degrevlex_benchmarks(ground)

#=
cyclic 7
  0.163683 seconds (91.60 k allocations: 55.431 MiB, 5.30% gc time, 31.38% compilation time) 
cyclic 8
  1.690214 seconds (267.60 k allocations: 292.095 MiB, 8.76% gc time, 0.23% compilation time)
katsura 9
  0.208915 seconds (74.23 k allocations: 53.478 MiB)
katsura 10
  1.492179 seconds (195.74 k allocations: 176.967 MiB, 11.31% gc time)
katsura 11
  7.061800 seconds (515.75 k allocations: 617.090 MiB, 1.46% gc time)
noon 7
  0.307783 seconds (209.98 k allocations: 77.458 MiB, 41.56% gc time)
noon 8
  1.665348 seconds (945.53 k allocations: 353.220 MiB, 14.79% gc time)
noon 9
 13.829567 seconds (4.30 M allocations: 2.044 GiB, 3.35% gc time)
eco 10
  0.115604 seconds (72.95 k allocations: 41.181 MiB)
eco 11
  0.538119 seconds (188.62 k allocations: 111.055 MiB)
eco 12
  2.762637 seconds (553.62 k allocations: 381.915 MiB)
=#
