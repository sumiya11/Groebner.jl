
import Singular
import AbstractAlgebra

using Base.Threads

import Pkg
using CpuId
using Printf
using BenchmarkTools

using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BENCHMARK_SAMPLES = 1

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

    systems = benchmark_systems(ground)

    resulting_md = ""
    resulting_md *= "|System|Singular|\n"

    for (name, system) in systems
        resulting_md *= "|$name|"
        run_system_singular(name, system, resulting_md)
        resulting_md *= "\n"
    end

    # resulting_md *= "\n*Groebner parameters:*\n"
    # resulting_md *= string(GROEBNER_PARAMS)
    resulting_md *= "\n*Benchmarking environment:*\n\n"
    resulting_md *= "* Total RAM (Mb): $(Sys.total_memory() / 2^20)\n"
    resulting_md *= "* Processor: $(cpubrand())\n"
    resulting_md *= "* Julia version: $(VERSION)\n\n"
    resulting_md *= "Versions of the dependencies:\n\n"

    open("singular_benchmark_result.md", "w") do io
        write(io, resulting_md)
    end
end

ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)
