
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
    ("cyclic 12", Groebner.cyclicn(7, ground=ground)),
    ("cyclic 13", Groebner.cyclicn(8, ground=ground)),
    ("cyclic 14", Groebner.cyclicn(9, ground=ground)),
    ("root 12", Groebner.rootn(12, ground=ground)),
    ("root 13", Groebner.rootn(13, ground=ground)),
    ("root 14", Groebner.rootn(14, ground=ground)),
        ("katsura 9", Groebner.katsuran(9, ground=ground)),
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
cyclic 12
67.366 ms (179903 allocations: 27.04 MiB)
cyclic 13
216.331 ms (451370 allocations: 67.29 MiB)
katsura 9
104.752 ms (55381 allocations: 17.84 MiB)
noon 6
31.177 ms (86994 allocations: 19.41 MiB)
noon 7
235.416 ms (419391 allocations: 89.92 MiB)
eco 10
182.549 ms (137522 allocations: 41.11 MiB)
ku 10
989.900 μs (5794 allocations: 2.52 MiB)
kinema
2.943 ms (6727 allocations: 3.27 MiB)

reduced
cyclic 12
397.781 ms (401132 allocations: 70.27 MiB)
cyclic 13
 1.465 s (1032878 allocations: 186.29 MiB)
katsura 9
 120.634 ms (56825 allocations: 18.90 MiB)
noon 6
 33.456 ms (89538 allocations: 20.19 MiB)
noon 7
 248.611 ms (426594 allocations: 93.01 MiB)
eco 10
 189.309 ms (139554 allocations: 42.71 MiB)
ku 10
 996.400 μs (5885 allocations: 2.53 MiB)
kinema
 2.715 ms (7193 allocations: 3.36 MiB)
=#
