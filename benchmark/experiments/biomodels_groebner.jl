
if !isdefined(Main, :Groebner)
    import Groebner
end

include("../systems/biomodels/parser.jl")

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1

function benchmark_system_my(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    # @btime gb = Groebner.groebner($system, linalg=:prob)
    stats = @timed @time Groebner.groebner(system, linalg=:prob, ordering=Groebner.DegRevLex())
    stats
end

function run_f4_biomodels_benchmarks(systems, threshold)
    println()
    goodnames = []
    for (name, system) in systems
        println("$name")
        if name in skip
            continue
        end
        stats = benchmark_system_my(system)
        if stats.time > threshold
            push!(goodnames, (name, stats.time))
        end
    end
    goodnames
end

skip = [
    "BIOMD0000000085",
    "BIOMD0000000086",
]

systems = read_BIOMDs(10:50)

goodnames = run_f4_biomodels_benchmarks(systems, 5.0)

#=

=#
