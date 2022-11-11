
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

    @btime Groebner.groebner($system, reduced=false, linalg=:prob)

end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
    ("cyclic 6", reverse(Groebner.cyclicn(6, ground=ground))),
    ("cyclic 7", reverse(Groebner.cyclicn(7, ground=ground))),
        ("katsura 7", reverse(Groebner.katsuran(7, ground=ground))),
        ("katsura 8", reverse(Groebner.katsuran(8, ground=ground))),
        ("katsura 9", reverse(Groebner.katsuran(9, ground=ground))),
        ("noon 6"    ,Groebner.noonn(6, ground=ground)),
        ("noon 7"    ,Groebner.noonn(7, ground=ground)),
        ("noon 8"    ,Groebner.noonn(8, ground=ground)),
        ("eco 10"    ,Groebner.eco10(ground=ground)),
        ("eco 11"    ,Groebner.eco11(ground=ground))
    ]

    for (name, system) in systems
        println("$name")
        benchmark_system_my(system)
    end
end

println()
ground = AbstractAlgebra.QQ
run_f4_ff_degrevlex_benchmarks(ground)


#=

qq:

cyclic 6
  16.083 ms (41950 allocations: 7.84 MiB)
cyclic 7
  4.348 s (3224917 allocations: 778.54 MiB)
katsura 7
  300.875 ms (395295 allocations: 67.01 MiB)
katsura 8
  2.269 s (2073004 allocations: 384.21 MiB)
katsura 9
  10.807 s (6798833 allocations: 1.07 GiB)
noon 6
  99.882 ms (293406 allocations: 42.26 MiB)
noon 7
  683.802 ms (1502368 allocations: 199.92 MiB)
noon 8
  7.004 s (9184113 allocations: 1.24 GiB)
eco 10
  651.303 ms (1018972 allocations: 208.98 MiB)
eco 11
  7.053 s (6136001 allocations: 1.22 GiB)

=#

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
noon 9
  12.392 s (4306578 allocations: 1.98 GiB)
eco 10
  66.734 ms (72768 allocations: 38.28 MiB)
eco 11
  320.164 ms (185723 allocations: 97.99 MiB)

=#
