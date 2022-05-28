
if !isdefined(Main, :Groebner)
    import Groebner
end

include("MQ/parser.jl")

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

    @btime gb = Groebner.groebner($system)
    # gb = Groebner.groebner(system, reduced=false)
end

function run_f4_MQ_benchmarks()
    systems = [
        ("mq_10_20_2_0", read_MQ_GF("mq_n10_m20_p2_s0")),
        ("mq_10_20_2_1", read_MQ_GF("mq_n10_m20_p2_s1"))
    ]

    println()
    for (name, system) in systems
        println("$name")
        # benchmark_system_my(system)
        benchmark_system_my(system)
    end
end

run_f4_MQ_benchmarks()

#=

=#
