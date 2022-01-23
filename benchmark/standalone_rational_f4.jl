
if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @btime gb = Groebner.groebner($system, reduced=false)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 11", Groebner.rootn(11, ground=ground)),
        ("cyclic 12", Groebner.rootn(12, ground=ground)),
        ("noon 5"    ,Groebner.noonn(5, ground=ground)),
        ("noon 6"    ,Groebner.noonn(6, ground=ground)),
    ]

    for (name, system) in systems
        println("$name")
        benchmark_system_my(system)
    end
end

ground = AbstractAlgebra.QQ
run_f4_ff_degrevlex_benchmarks(ground)

#=
cyclic 11
51.978 ms (218960 allocations: 21.78 MiB)
cyclic 12
168.396 ms (506977 allocations: 53.18 MiB)
noon 5
23.182 ms (129254 allocations: 13.04 MiB)
noon 6
227.689 ms (830899 allocations: 88.88 MiB)
=#
