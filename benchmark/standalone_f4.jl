
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
ground = AbstractAlgebra.GF(2)
run_f4_ff_degrevlex_benchmarks(ground)

#=
modulo 2^31 - 1

cyclic 7
  265.497 ms (96747 allocations: 54.72 MiB)
cyclic 8
  5.506 s (435423 allocations: 276.82 MiB)
root 12
  69.753 ms (179996 allocations: 31.06 MiB)
root 13
  235.030 ms (451468 allocations: 63.35 MiB)
katsura 9
  784.184 ms (123589 allocations: 47.47 MiB)
noon 7
  220.374 ms (335526 allocations: 75.72 MiB)
noon 8
  2.173 s (1501614 allocations: 378.00 MiB)
eco 10
  197.265 ms
=#
