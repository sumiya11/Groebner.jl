
if !isdefined(Main, :Groebner)
    using Groebner
end

import AbstractAlgebra, Primes
using BenchmarkTools
using Logging

Groebner.invariants_enabled() = false

# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 6
BenchmarkTools.DEFAULT_PARAMETERS.samples = 3

Groebner.logging_enabled() = false
Groebner.invariants_enabled() = true

function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    @btime Groebner.groebner($system)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 7", reverse(Groebner.cyclicn(7, k=ground))),
        ("cyclic 8", reverse(Groebner.cyclicn(8, k=ground))),
        ("katsura 9", reverse(Groebner.katsuran(9, k=ground))),
        ("katsura 10", reverse(Groebner.katsuran(10, k=ground))),
        ("katsura 11", reverse(Groebner.katsuran(11, k=ground))),
        ("noon 7", Groebner.noonn(7, k=ground)),
        ("noon 8", Groebner.noonn(8, k=ground)),
        # ("noon 9", Groebner.noonn(9, k=ground)),
        ("eco 10", Groebner.eco10(k=ground)),
        ("eco 11", Groebner.eco11(k=ground))
    ]

    for (name, system) in systems
        println("$name")
        benchmark_system_my(system)
    end
end

p1 = 2^31 - 1
p2 = 2^29 + 11
p3 = 2^27 + 29
s = Groebner.katsuran(11, ordering=:degrevlex, k=AbstractAlgebra.GF(p3))

@profview gb1 = Groebner.groebner(s);

println()
ground = AbstractAlgebra.GF(1048583)
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
