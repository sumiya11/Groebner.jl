
if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10
# BenchmarkTools.DEFAULT_PARAMETERS.samples = 2


function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # println("length = $(length(system)), nvars = $(AbstractAlgebra.nvars(AbstractAlgebra.parent(system[1])))")
    # println("maxdeg = $(maximum(AbstractAlgebra.total_degree ∘ AbstractAlgebra.leading_monomial, system))")

    @btime Groebner.groebner($system, reduced=true)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
    ("cyclic 10", Groebner.rootn(10, ground=ground)),
    ("cyclic 11", Groebner.rootn(11, ground=ground)),
    ("cyclic 12", Groebner.rootn(12, ground=ground)),
    ("cyclic 13", Groebner.rootn(13, ground=ground)),
        ("katsura 6", Groebner.katsura6(ground=ground)),
        ("katsura 9", Groebner.katsura9(ground=ground)),
        ("katsura 10", Groebner.katsura10(ground=ground)),
        ("noon 4"    ,Groebner.noonn(4, ground=ground)),
        ("noon 5"    ,Groebner.noonn(5, ground=ground)),
        ("noon 6"    ,Groebner.noonn(6, ground=ground)),
        ("noon 7"    ,Groebner.noonn(7, ground=ground)),
        ("eco 7"    ,Groebner.eco7(ground=ground)),
        ("eco 10"    ,Groebner.eco10(ground=ground)),
        ("eco 11"    ,Groebner.eco11(ground=ground)),
        ("ku 10"    ,Groebner.ku10(ground=ground)),
        ("s9_1"    ,Groebner.s9_1(ground=ground)),
        ("kinema"    ,Groebner.kinema(ground=ground)),
        ("ojika4_d1R2_d2R5"    ,Groebner.ojika4_d1R2_d2R5(ground=ground))
    ]

    println()

    for (name, system) in systems
        println("$name")
        benchmark_system_my(system)
    end
end

function big_run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        #("cyclic 14", Groebner.rootn(14, ground=ground)),
        #("cyclic 15", Groebner.rootn(15, ground=ground)),
        #("katsura 11", Groebner.katsura11(ground=ground)),
        #("katsura 12", Groebner.katsura12(ground=ground)),
        #("katsura 13", Groebner.katsura13(ground=ground)),
        ("noon 7"    ,Groebner.noonn(7, ground=ground)),
        ("noon 8"    ,Groebner.noonn(8, ground=ground)),
        ("eco 11"    ,Groebner.eco11(ground=ground)),
        ("eco 12"    ,Groebner.eco12(ground=ground)),
        #("eco 13"    ,Groebner.eco13(ground=ground)),
    ]

    println()

    for (name, system) in systems
        println("$name")
        benchmark_system_my(system, tablesize)
    end
end

println()
ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)

#=
2^6:
cyclic 10
  25.030 ms (61954 allocations: 10.70 MiB)
cyclic 11
  84.327 ms (156817 allocations: 25.69 MiB)
cyclic 12
  352.658 ms (401131 allocations: 70.27 MiB)
cyclic 13
  1.492 s (1032877 allocations: 186.29 MiB)
katsura 6
  3.615 ms (8778 allocations: 2.07 MiB)
katsura 9
  112.370 ms (56824 allocations: 18.90 MiB)
katsura 10
  765.672 ms (152269 allocations: 62.73 MiB)
noon 4
  1.133 ms (4518 allocations: 831.12 KiB)
noon 5
  5.000 ms (19465 allocations: 3.80 MiB)
noon 6
  36.483 ms (89537 allocations: 20.19 MiB)
noon 7
  253.887 ms (426593 allocations: 93.01 MiB)
eco 7
  2.122 ms (8617 allocations: 2.13 MiB)
eco 10
  196.016 ms (139553 allocations: 42.71 MiB)
eco 11
  1.302 s (394404 allocations: 136.15 MiB)
ku 10
  1.146 ms (5884 allocations: 2.53 MiB)
s9_1
  495.600 μs (1901 allocations: 1.82 MiB)
kinema
  2.961 ms (7192 allocations: 3.36 MiB)
ojika4_d1R2_d2R5
  311.500 μs (2154 allocations: 354.58 KiB)

2^4:

2^7:

=#
