
if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 10


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
    ("cyclic 12", Groebner.rootn(12, ground=ground)),
    ("cyclic 13", Groebner.rootn(13, ground=ground)),
        ("katsura 9", Groebner.katsura9(ground=ground)),
        ("noon 6"    ,Groebner.noonn(6, ground=ground)),
        ("noon 7"    ,Groebner.noonn(7, ground=ground)),
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
66.675 ms (179922 allocations: 27.62 MiB)
cyclic 13
221.970 ms (451391 allocations: 68.78 MiB)
katsura 9
319.554 ms (55390 allocations: 17.93 MiB)
noon 6
34.842 ms (87010 allocations: 20.69 MiB)
noon 7
265.669 ms (419413 allocations: 91.78 MiB)
eco 10
470.273 ms (137548 allocations: 41.31 MiB)
=#
