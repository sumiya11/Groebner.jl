
import Pkg
Pkg.add("Groebner")
Pkg.add("AbstractAlgebra")
Pkg.add("BenchmarkTools")

using Groebner
import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    gb = Groebner.groebner(system)

    @btime gb = Groebner.groebner($system, reduced=true)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("reimer 6", Groebner.reimern(6, ground=ground)),
        # ("reimer 7", Groebner.reimern(7, ground=ground)),
        # ("reimer 8", Groebner.reimern(8, ground=ground)),
        # ("cyclic 12", Groebner.rootn(12, ground=ground)),
        # ("cyclic 13", Groebner.rootn(13, ground=ground)),
        # ("katsura 9",Groebner.katsura9(ground=ground)),
        # ("katsura 10",Groebner.katsura10(ground=ground)),
        # ("eco 10",Groebner.eco10(ground=ground)),
        # ("eco 11",Groebner.eco11(ground=ground)),
        # ("noon 7"    ,Groebner.noonn(7, ground=ground)),
        # ("noon 8"    ,Groebner.noonn(8, ground=ground))
    ]

    for (name, system) in systems
        println("$name")
        println("my")
        benchmark_system_my(system)
    end
end

ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)
