
include("../src/Groebner.jl")


if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 6
BenchmarkTools.DEFAULT_PARAMETERS.samples = 5


function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @btime Groebner.groebner($system, reduced=false)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
    ("cyclic 7", Groebner.cyclicn(7, ground=ground)),
    ("cyclic 8", Groebner.cyclicn(8, ground=ground)),
    ("root 13", Groebner.rootn(13, ground=ground)),
    ("root 14", Groebner.rootn(14, ground=ground)),
    ("root 15", Groebner.rootn(15, ground=ground)),
        ("katsura 9", Groebner.katsuran(9, ground=ground)),
        ("katsura 10", Groebner.katsuran(10, ground=ground)),
        ("katsura 11", Groebner.katsuran(11, ground=ground)),
        ("noon 7"    ,Groebner.noonn(7, ground=ground)),
        ("noon 8"    ,Groebner.noonn(8, ground=ground)),
        ("noon 9"    ,Groebner.noonn(9, ground=ground)),
        ("eco 10"    ,Groebner.eco10(ground=ground)),
        ("eco 11"    ,Groebner.eco11(ground=ground)),
        ("eco 12"    ,Groebner.eco12(ground=ground))
    ]

    for (name, system) in systems
        println("$name")
        benchmark_system_my(system)
    end
end

println()
ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)

#=
cyclic 7
212.187 ms (229506 allocations: 79.33 MiB)
cyclic 8
4.631 s (1362066 allocations: 455.37 MiB)
root 12
73.354 ms (179930 allocations: 32.39 MiB)
root 13
248.371 ms (451378 allocations: 67.86 MiB)
root 14
726.063 ms (1144701 allocations: 176.53 MiB)
katsura 9
882.256 ms (149540 allocations: 59.57 MiB)
noon 7
249.928 ms (419411 allocations: 91.51 MiB)
noon 8
1.681 s (2074539 allocations: 492.01 MiB)
eco 10
179.106 ms (137529 allocations: 41.41 MiB)
=#
