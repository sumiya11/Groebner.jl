
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
  216.897 ms (96622 allocations: 70.60 MiB)
cyclic 8
  4.521 s (435192 allocations: 380.99 MiB)
root 12
  74.583 ms (179896 allocations: 32.40 MiB)
root 13
  242.372 ms (451340 allocations: 67.87 MiB)
katsura 9
  761.210 ms (123510 allocations: 59.33 MiB)
noon 7
  248.028 ms (335407 allocations: 86.89 MiB)
noon 8
  1.960 s (1501606 allocations: 444.08 MiB)
eco 10
  191.307 ms (119125 allocations: 41.66 MiB)
=#
