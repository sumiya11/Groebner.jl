
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
Groebner.invariants_enabled() = false

function benchmark_system_my(system)
    Groebner.groebner([system[1]])

    @btime Groebner.groebner($system)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 7", reverse(Groebner.cyclicn(7, k=ground, internal_ordering=:degrevlex))),
        ("cyclic 8", reverse(Groebner.cyclicn(8, k=ground, internal_ordering=:degrevlex))),
        (
            "katsura 9",
            reverse(Groebner.katsuran(9, k=ground, internal_ordering=:degrevlex))
        ),
        (
            "katsura 10",
            reverse(Groebner.katsuran(10, k=ground, internal_ordering=:degrevlex))
        ),
        (
            "katsura 11",
            reverse(Groebner.katsuran(11, k=ground, internal_ordering=:degrevlex))
        ),
        ("noon 7", Groebner.noonn(7, k=ground, internal_ordering=:degrevlex)),
        ("noon 8", Groebner.noonn(8, k=ground, internal_ordering=:degrevlex)),
        # ("noon 9", Groebner.noonn(9, k=ground, internal_ordering=:degrevlex)),
        ("eco 10", Groebner.eco10(k=ground, internal_ordering=:degrevlex)),
        ("eco 11", Groebner.eco11(k=ground, internal_ordering=:degrevlex))
    ]

    for (name, system) in systems
        println("$name")
        benchmark_system_my(system)
    end
end

p1 = 2^30 + 3
s = Groebner.katsuran(11, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p1))

# @profview gb1 = Groebner.groebner(s);

println()
ground = AbstractAlgebra.GF(2^30 + 3)
run_f4_ff_degrevlex_benchmarks(ground)

#=

gf:

cyclic 7
  80.674 ms (32571 allocations: 39.23 MiB)
cyclic 8
  1.078 s (158547 allocations: 162.62 MiB)
katsura 9
  126.103 ms (34041 allocations: 43.78 MiB)
katsura 10
  708.583 ms (86954 allocations: 147.93 MiB)
katsura 11
  4.486 s (217870 allocations: 518.37 MiB)
noon 7
  127.639 ms (92085 allocations: 51.93 MiB)
noon 8
  1.097 s (429005 allocations: 283.11 MiB)
eco 10
  56.632 ms (37665 allocations: 29.22 MiB)
eco 11
  334.697 ms (125935 allocations: 83.94 MiB)

=#
