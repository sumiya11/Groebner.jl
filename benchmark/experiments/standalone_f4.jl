
include("../../src/Groebner.jl")

if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 6
BenchmarkTools.DEFAULT_PARAMETERS.samples = 3

Groebner.logging_enabled() = false
Groebner.invariants_enabled() = true

function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @time Groebner.groebner(system)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 7", reverse(Groebner.cyclicn(7, ground=ground))),
        ("cyclic 8", reverse(Groebner.cyclicn(8, ground=ground))),
        ("katsura 9", reverse(Groebner.katsuran(9, ground=ground))),
        ("katsura 10", reverse(Groebner.katsuran(10, ground=ground))),
        ("katsura 11", reverse(Groebner.katsuran(11, ground=ground))),
        ("noon 7", Groebner.noonn(7, ground=ground)),
        ("noon 8", Groebner.noonn(8, ground=ground)),
        # ("noon 9", Groebner.noonn(9, ground=ground)),
        ("eco 10", Groebner.eco10(ground=ground)),
        ("eco 11", Groebner.eco11(ground=ground))
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

gf:

cyclic 7
  95.218 ms (58110 allocations: 51.84 MiB)
cyclic 8
  1.415 s (268517 allocations: 278.05 MiB)
katsura 9
  155.376 ms (75093 allocations: 50.72 MiB)
katsura 10
  943.326 ms (204033 allocations: 173.79 MiB)
katsura 11
  6.671 s (525176 allocations: 612.59 MiB)
noon 7
  194.497 ms (209483 allocations: 74.25 MiB)
noon 8
  1.327 s (950722 allocations: 378.90 MiB)
eco 10
  66.734 ms (72768 allocations: 38.28 MiB)
eco 11
  320.164 ms (185723 allocations: 97.99 MiB)

=#
