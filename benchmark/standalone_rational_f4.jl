
if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BenchmarkTools.DEFAULT_PARAMETERS.samples = 5

function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @btime gb = Groebner.groebner($system, reduced=false)
    # gb = Groebner.groebner(system, reduced=false)
end

function benchmark_system_my_certify(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @btime gb = Groebner.groebner($system, reduced=false, linalg=:prob)
    # gb = Groebner.groebner(system, reduced=false)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("siwr", read_SIWR()),
        ("sear", read_SEAIJRC()),
        ("root 11", Groebner.rootn(11, ground=ground)),
        ("noon 7"    ,Groebner.noonn(7, ground=ground)),
        ("eco 10"    ,Groebner.eco10(ground=ground)),
        ("katsura 8"    ,Groebner.katsuran(8, ground=ground))
    ]

    println()
    for (name, system) in systems
        println("$name")
        # benchmark_system_my(system)
        benchmark_system_my_certify(system)
    end
end

ground = AbstractAlgebra.QQ
run_f4_ff_degrevlex_benchmarks(ground)

#=
root 11
  55.078 ms (183751 allocations: 23.39 MiB)
noon 7
  805.825 ms (1773797 allocations: 276.06 MiB)
eco 10
  911.046 ms (1324278 allocations: 294.65 MiB)
katsura 8
  3.124 s (1794620 allocations: 312.78 MiB)
=#
