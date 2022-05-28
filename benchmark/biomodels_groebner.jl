
if !isdefined(Main, :Groebner)
    import Groebner
end

include("biomodels/parser.jl")

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BenchmarkTools.DEFAULT_PARAMETERS.samples = 3

function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    @btime gb = Groebner.groebner($system)
    # gb = Groebner.groebner(system, reduced=false)
end

function run_f4_biomodels_benchmarks()
    systems = read_BIOMDs(5:20)

    println()
    for (name, system) in systems
        println("$name")
        # benchmark_system_my(system)
        benchmark_system_my(system)
    end
end

run_f4_biomodels_benchmarks()

#=

=#
