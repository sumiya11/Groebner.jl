
import Singular
import AbstractAlgebra

using Base.Threads

import Pkg
using CpuId
using Printf
using BenchmarkTools
include("../../src/Groebner.jl")

using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BENCHMARK_SAMPLES = 1
GROEBNER_PARAMS = (linalg=:prob, )

function benchmark_system_groebner(system)
    system = Groebner.change_ordering(system, :degrevlex)
    # Eliminate compilation times
    global GROEBNER_PARAMS
    global BENCHMARK_SAMPLES
    println("Basis size = ", length(Groebner.groebner(system; GROEBNER_PARAMS...)))
    # run benchmarks
    bench = @benchmarkable Groebner.groebner($system; $GROEBNER_PARAMS...) samples=BENCHMARK_SAMPLES
    median(run(bench))
end

function benchmark_system_singular(system)
    R = AbstractAlgebra.parent(system[1])
    modulo = AbstractAlgebra.characteristic(R)
    n = AbstractAlgebra.nvars(R)
    ground_s = Singular.N_ZpField(modulo)
    R_s, _ = Singular.PolynomialRing(ground_s, ["x$i" for i in 1:n], ordering=:degrevlex)
    system_s = map(
        f -> AbstractAlgebra.change_base_ring(
                    ground_s,
                    AbstractAlgebra.map_coefficients(c -> ground_s(c.d), f),
                    parent=R_s),
        system)
    ideal_s = Singular.Ideal(R_s, system_s)
    # Eliminate compilation times
    Singular.std(Singular.Ideal(R_s, [system_s[1]]))
    global BENCHMARK_SAMPLES
    bench = @benchmarkable  Singular.std($ideal_s) samples=BENCHMARK_SAMPLES
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

function run_system_singular(name, system, resulting_md)
    x = benchmark_system_singular(system)
    println("---------------\n$name")
    println("singular")
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
    systems = [
        # ("cyclic 7",  Groebner.cyclicn(7, ground=ground)), 
        # ("cyclic 8",  Groebner.cyclicn(8, ground=ground)), 
        # ("cyclic 9",  Groebner.cyclicn(9, ground=ground)), 
        ("katsura 10",Groebner.katsuran(10, ground=ground)), 
        ("katsura 11",Groebner.katsuran(11, ground=ground)), 
        ("katsura 12",Groebner.katsuran(12, ground=ground)), 
        ("eco 11",    Groebner.eco11(ground=ground)),         
        ("eco 12",    Groebner.eco12(ground=ground)),         
        ("eco 13",    Groebner.eco13(ground=ground)),         
        ("noon 7",    Groebner.noonn(7, ground=ground)),  
        ("noon 8",    Groebner.noonn(8, ground=ground)),  
        ("noon 9",    Groebner.noonn(9, ground=ground)),   
        ("henrion 5", Groebner.henrion5(ground=ground)),  
        ("henrion 6", Groebner.henrion6(ground=ground)),  
        ("henrion 7", Groebner.henrion7(ground=ground)),  
        ("reimer 6",  Groebner.reimern(6, ground=ground)), 
        ("reimer 7",  Groebner.reimern(7, ground=ground)), 
        ("reimer 8",  Groebner.reimern(8, ground=ground)), 
    ]

    resulting_md = ""
    resulting_md *= "|System|Groebner.jl|Singular|\n"

    for (name, system) in systems
        resulting_md *= "|$name|"
        run_system_groebner(name, system, resulting_md)
        # run_system_singular(name, system, resulting_md)
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
        if deps[uid].version != nothing
            resulting_md *= "* $s : $(deps[uid].version)\n"
        end
    end

    open("benchmark_result.md", "w") do io
        write(io, resulting_md)
    end
end

ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)
