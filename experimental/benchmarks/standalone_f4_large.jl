
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
        ("cyclic 7", Groebner.cyclicn(7, k=ground)),
        ("cyclic 8", Groebner.cyclicn(8, k=ground)),
        ("cyclic 9", Groebner.cyclicn(9, k=ground)),
        ("katsura 9", Groebner.katsuran(9, k=ground)),
        ("katsura 10", Groebner.katsuran(10, k=ground)),
        ("katsura 11", Groebner.katsuran(11, k=ground)),
        ("katsura 12", Groebner.katsuran(12, k=ground)),
        ("noon 7", Groebner.noonn(7, k=ground)),
        ("noon 8", Groebner.noonn(8, k=ground)),
        ("noon 9", Groebner.noonn(9, k=ground)),
        ("eco 11", Groebner.eco11(k=ground)),
        ("eco 12", Groebner.eco12(k=ground)),
        ("eco 12", Groebner.eco13(k=ground))
    ]

    for (name, system) in systems
        print("$name -- ")
        t = benchmark_system_my(system)
        println("$t s")
    end
end

ground = AbstractAlgebra.GF(2^30 + 3)
run_f4_ff_degrevlex_benchmarks(ground)
