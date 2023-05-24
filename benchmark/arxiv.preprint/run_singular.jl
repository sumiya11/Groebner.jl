
# Load all required packages
import Pkg
Pkg.add("Groebner")
Pkg.add("AbstractAlgebra")
Pkg.add("BenchmarkTools")
Pkg.add("Singular")

# Import the packages
using Groebner
import Singular
import AbstractAlgebra
using BenchmarkTools
using Logging
include((@__DIR__) * "/generate/benchmark_systems.jl")
global_logger(ConsoleLogger(stderr, Logging.Error))

# Benchmark the given system
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
            parent=R_s
        ),
        system
    )
    ideal_s = Singular.Ideal(R_s, system_s)
    # compile:
    Singular.std(Singular.Ideal(R_s, [system_s[1]]))
    # run the actual benchmark
    @btime gb = Singular.std($ideal_s, complete_reduction=true)
end

function run_f4_ff_degrevlex_benchmarks(flag)
    if flag
        ground = AbstractAlgebra.GF(2^31 - 1)
        systems = benchmark_systems_ff(ground)
    else
        ground = AbstractAlgebra.QQ
        systems = benchmark_systems_qq(ground)
    end

    println("Running Singular benchmarks over $ground")
    for (name, system) in systems
        println("==================")
        println("System $name:")
        benchmark_system_singular(system)
    end
end

run_f4_ff_degrevlex_benchmarks(true)
run_f4_ff_degrevlex_benchmarks(false)
