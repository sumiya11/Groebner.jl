
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
sqrt(nlow / 3)
cyclic 7 -- 0.070263489 s
cyclic 8 -- 1.14754937 s
cyclic 9 -- 95.630521758 s
katsura 9 -- 0.121808741 s
katsura 10 -- 0.647219826 s
katsura 11 -- 4.292973661 s
katsura 12 -- 31.899829963 s
noon 7 -- 0.140998391 s
noon 8 -- 1.506614076 s
noon 9 -- 15.763489396 s
eco 11 -- 0.308983478 s
eco 12 -- 1.780714725 s
eco 12 -- 8.237998759 s
=#
