
if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10
BenchmarkTools.DEFAULT_PARAMETERS.samples = 2


function benchmark_system_my(system, tablesize)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]], tablesize)

    println("length = $(length(system)), nvars = $(AbstractAlgebra.nvars(AbstractAlgebra.parent(system[1])))")
    println("maxdeg = $(maximum(AbstractAlgebra.total_degree ∘ AbstractAlgebra.leading_monomial, system))")

    @btime Groebner.groebner($system, $tablesize, reduced=true)
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
        for tablesize in (2^13, 2^14, 2^15, 2^16)
            println("$name - $tablesize")
            benchmark_system_my(system, tablesize)
        end
        println("########\n")
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
        for tablesize in (2^12, 2^14, 2^16)
            println("$name - $tablesize")
            benchmark_system_my(system, tablesize)
        end
        println("########\n")
    end
end

println()
ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)

#=

=#
