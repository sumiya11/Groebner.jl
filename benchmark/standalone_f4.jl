
if !isdefined(Main, :FastGroebner)
    import FastGroebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

function benchmark_system_my(system)
    system = FastGroebner.change_ordering(system, :degrevlex)
    FastGroebner.groebner([system[1]])
    @btime gb = FastGroebner.groebner($system, reduced=false)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 12", FastGroebner.rootn(12, ground=ground)),
        ("cyclic 13", FastGroebner.rootn(13, ground=ground)),
        ("katsura 9", FastGroebner.katsura9(ground=ground)),
        ("katsura 10",FastGroebner.katsura10(ground=ground)),
        ("noon 6"    ,FastGroebner.noonn(6, ground=ground)),
        ("noon 7"    ,FastGroebner.noonn(7, ground=ground))
    ]

    for (name, system) in systems
        println("$name")
        benchmark_system_my(system)
    end
end

ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)

#=
cyclic 12
  78.942 ms (184266 allocations: 29.24 MiB)
cyclic 13
  272.030 ms (459882 allocations: 74.68 MiB)
katsura 9
  370.388 ms (82529 allocations: 22.35 MiB)
katsura 10
  2.969 s (250068 allocations: 74.80 MiB)
noon 6
  42.705 ms (101234 allocations: 21.19 MiB)
noon 7
  336.392 ms (494526 allocations: 100.91 MiB)
=#
