
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

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 7", Groebner.cyclicn(7, ground=ground)),
        ("cyclic 8", Groebner.cyclicn(8, ground=ground)),
        ("katsura 10", Groebner.katsuran(10, ground=ground)),
        ("katsura 11", Groebner.katsuran(11, ground=ground)),
        ("eco 12", Groebner.eco12(ground=ground)),
        ("eco 13", Groebner.eco13(ground=ground)),
        ("noon 7", Groebner.noonn(7, ground=ground)),
        ("noon 8", Groebner.noonn(8, ground=ground))
    ]

    println("Running Singular benchmarks")
    for (name, system) in systems
        println("==================")
        println("System $name:")
        benchmark_system_singular(system)
    end
end

ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)
