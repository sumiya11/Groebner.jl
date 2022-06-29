
include("standard/parser.jl")

import AbstractAlgebra
using BenchmarkTools
using Logging

global_logger(ConsoleLogger(stderr, Logging.Error))

BenchmarkTools.DEFAULT_PARAMETERS.samples = 2

function benchmark_system(system)
    @btime Groebner.groebner($s, ordering=:degrevlex)
end

function run_for_gleb()
    fns = [
        "gbSIWR",
        "gbSEAIJRC"
    ]

    for fn in fns
        system = load_system(fn)
        benchmark_system(system)
    end
end
