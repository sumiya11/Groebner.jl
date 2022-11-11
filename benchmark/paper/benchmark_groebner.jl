
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

BENCHMARK_SAMPLES = 2
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
    # stid_info = deps[findfirst(x -> x.name == "Groebner", deps)]
    # for (s, uid) in stid_info.dependencies
    #     if deps[uid].version !== nothing
    #         resulting_md *= "* $s : $(deps[uid].version)\n"
    #     end
    # end

    open("groebner_benchmark_result.md", "w") do io
        write(io, resulting_md)
    end
end

ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)

# with prob
#=
Basis size = 209
---------------
cyclic 7
groebner
TrialEstimate(116.678 ms)
---------------
Basis size = 372
---------------
cyclic 8
groebner
TrialEstimate(2.200 s)
---------------
Basis size = 1344
---------------
cyclic 9
groebner
TrialEstimate(253.007 s)
---------------
Basis size = 537
---------------
katsura 10
groebner
TrialEstimate(1.072 s)
---------------
Basis size = 1050
---------------
katsura 11
groebner
TrialEstimate(7.890 s)
---------------
Basis size = 2091
---------------
katsura 12
groebner
TrialEstimate(55.039 s)
---------------
Basis size = 389
---------------
eco 11
groebner
TrialEstimate(365.310 ms)
---------------
Basis size = 743
---------------
eco 12
groebner
TrialEstimate(2.294 s)
---------------
Basis size = 1465
---------------
eco 13
groebner
TrialEstimate(11.412 s)
---------------
Basis size = 495
---------------
noon 7
groebner
TrialEstimate(144.760 ms)
---------------
Basis size = 1338
---------------
noon 8
groebner
TrialEstimate(1.266 s)
---------------
Basis size = 3682
---------------
noon 9
groebner
TrialEstimate(11.588 s)
---------------
Basis size = 33
---------------
henrion 5
groebner
TrialEstimate(2.105 ms)
---------------
Basis size = 90
---------------
henrion 6
groebner
TrialEstimate(32.866 ms)
---------------
Basis size = 415
---------------
henrion 7
groebner
TrialEstimate(3.205 s)
---------------
Basis size = 95
---------------
reimer 6
groebner
TrialEstimate(47.753 ms)
---------------
Basis size = 227
---------------
reimer 7
groebner
TrialEstimate(928.748 ms)
---------------
Basis size = 612
---------------
reimer 8
groebner
TrialEstimate(26.609 s)
---------------
542
=#

# without prob
#=
Basis size = 209
---------------
cyclic 7
groebner
TrialEstimate(165.265 ms)
---------------
Basis size = 372
---------------
cyclic 8
groebner
TrialEstimate(3.783 s)
---------------
Basis size = 1344
---------------
cyclic 9
groebner
TrialEstimate(531.727 s)
---------------
Basis size = 537
---------------
katsura 10
groebner
TrialEstimate(5.127 s)
---------------
Basis size = 1050
---------------
katsura 11
groebner
TrialEstimate(50.432 s)
---------------
Basis size = 2091
---------------
katsura 12
groebner
TrialEstimate(425.632 s)
---------------
Basis size = 389
---------------
eco 11
groebner
TrialEstimate(1.029 s)
---------------
Basis size = 743
---------------
eco 12
groebner
TrialEstimate(8.218 s)
---------------
Basis size = 1465
---------------
eco 13
groebner
TrialEstimate(76.818 s)
---------------
Basis size = 495
---------------
noon 7
groebner
TrialEstimate(143.092 ms)
---------------
Basis size = 1338
---------------
noon 8
groebner
TrialEstimate(1.225 s)
---------------
Basis size = 3682
---------------
noon 9
groebner
TrialEstimate(8.974 s)
---------------
Basis size = 33
---------------
henrion 5
groebner
TrialEstimate(2.376 ms)
---------------
Basis size = 90
---------------
henrion 6
groebner
TrialEstimate(40.457 ms)
---------------
Basis size = 415
---------------
henrion 7
groebner
TrialEstimate(5.423 s)
---------------
Basis size = 95
---------------
reimer 6
groebner
TrialEstimate(63.541 ms)
---------------
Basis size = 227
---------------
reimer 7
groebner
TrialEstimate(1.490 s)
---------------
Basis size = 612
---------------
reimer 8
groebner
TrialEstimate(52.370 s)
---------------
543
=#