
if !isdefined(Main, :Groebner)
    import Groebner
end

include("biomodels/parser.jl")

import Singular
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

function benchmark_system_singular(system)
    R = AbstractAlgebra.parent(system[1])
    n = AbstractAlgebra.nvars(R)
    ground_s = Singular.QQ
    R_s, _ = Singular.polynomial_ring(
        ground_s,
        ["x$i" for i in 1:n],
        internal_ordering=:degrevlex
    )

    system_s = map(
        f -> AbstractAlgebra.change_base_ring(
            ground_s,
            AbstractAlgebra.map_coefficients(
                c -> ground_s(numerator(c), denominator(c)),
                f
            ),
            parent=R_s
        ),
        system
    )

    ideal_s = Singular.Ideal(R_s, system_s)

    Singular.std(Singular.Ideal(R_s, [system_s[1]]))

    @btime Singular.std($ideal_s)
end

function run_biomodels_benchmarks()
    systems = read_BIOMDs(5:20)

    println()
    for (name, system) in systems
        println("$name")
        println("my")
        benchmark_system_my(system)
        println("singular")
        benchmark_system_singular(system)
    end
end

run_biomodels_benchmarks()

#=

=#
