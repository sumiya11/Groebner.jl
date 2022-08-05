
include("../../src/Groebner.jl")

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
  220.913 ms (96743 allocations: 54.72 MiB)
cyclic 8
  4.700 s (435419 allocations: 276.82 MiB)
root 12
  72.045 ms (179992 allocations: 31.06 MiB)
root 13
  229.356 ms (451464 allocations: 63.35 MiB)
katsura 9
  799.492 ms (123585 allocations: 47.47 MiB)
noon 7
  225.302 ms (335522 allocations: 75.72 MiB)
noon 8
  1.785 s (1501610 allocations: 378.00 MiB)
eco 10
  196.429 ms (119221 allocations: 33.93 MiB)
=#
