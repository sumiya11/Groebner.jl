
import AbstractAlgebra

using Base.Threads

import Pkg
using CpuId
using Printf
using BenchmarkTools
include("../../src/Groebner.jl")

include("benchmark_systems.jl")

using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BENCHMARK_SAMPLES = 1
GROEBNER_PARAMS = (linalg=:prob, )

function benchmark_system_groebner(system)
    system = Groebner.change_ordering(system, :degrevlex)
    # Reduce compilation times
    global GROEBNER_PARAMS
    global BENCHMARK_SAMPLES
    println("Basis size = ", length(Groebner.groebner(system; GROEBNER_PARAMS...)))
    # run benchmarks
    bench = @benchmarkable Groebner.groebner($system; $GROEBNER_PARAMS...) samples=BENCHMARK_SAMPLES
    median(run(bench))
end

function run_system_groebner(name, system, resulting_md)
    x = benchmark_system_groebner(system)
    println("---------------\n$name")
    println("groebner")
    println(x)
    println("---------------")
    resulting_md *= string(x) * "|"
end

#= 
    A part of this benchmarking script is taken 
    from StructuralIdentifiability.jl benchmarks
    https://github.com/SciML/StructuralIdentifiability.jl/blob/master/benchmarking/run_benchmarks.jl
=#

function run_f4_ff_degrevlex_benchmarks(ground)

    systems = benchmark_systems(ground)

    resulting_md = ""
    resulting_md *= "|System|Groebner.jl|\n"

    for (name, system) in systems
        resulting_md *= "|$name|"
        run_system_groebner(name, system, resulting_md)
        resulting_md *= "\n"
    end

    global GROEBNER_PARAMS
    resulting_md *= "\n*Groebner parameters:*\n"
    resulting_md *= string(GROEBNER_PARAMS)
    resulting_md *= "\n*Benchmarking environment:*\n\n"
    resulting_md *= "* Total RAM (Mb): $(Sys.total_memory() / 2^20)\n"
    resulting_md *= "* Processor: $(cpubrand())\n"
    resulting_md *= "* Julia version: $(VERSION)\n\n"
    resulting_md *= "Versions of the dependencies:\n\n"

    deps = Pkg.dependencies()
    stid_info = deps[findfirst(x -> x.name == "Groebner", deps)]
    for (s, uid) in stid_info.dependencies
        if deps[uid].version !== nothing
            resulting_md *= "* $s : $(deps[uid].version)\n"
        end
    end

    open("groebner_benchmark_result.md", "w") do io
        write(io, resulting_md)
    end
end

ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)
