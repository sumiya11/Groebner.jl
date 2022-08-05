
include("../src/Groebner.jl")


if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 6
BenchmarkTools.DEFAULT_PARAMETERS.samples = 4


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
    ("cyclic 7", reverse(Groebner.cyclicn(7, ground=ground))),
    ("cyclic 8", reverse(Groebner.cyclicn(8, ground=ground))),
    ("root 12", reverse(Groebner.rootn(12, ground=ground))),
    ("root 13", reverse(Groebner.rootn(13, ground=ground))),
    #("root 14", Groebner.rootn(14, ground=ground)),
        ("katsura 9", reverse(Groebner.katsuran(9, ground=ground))),
        ("noon 7"    ,Groebner.noonn(7, ground=ground)),
        ("noon 8"    ,Groebner.noonn(8, ground=ground)),
        ("eco 10"    ,Groebner.eco10(ground=ground))
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
  209.316 ms (229506 allocations: 79.33 MiB)
cyclic 8
  4.351 s (1362066 allocations: 455.37 MiB)
root 12
  61.705 ms (179930 allocations: 32.39 MiB)
root 13
  211.091 ms (451378 allocations: 67.86 MiB)
root 14
  809.867 ms (1144701 allocations: 176.53 MiB)
katsura 9
  801.246 ms (149540 allocations: 59.57 MiB)
noon 7
  214.115 ms (419411 allocations: 91.51 MiB)
noon 8
  1.842 s (2074539 allocations: 492.01 MiB)
eco 10
  179.854 ms (137529 allocations: 41.41 MiB)
=#
