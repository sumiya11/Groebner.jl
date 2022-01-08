
import Singular
import AbstractAlgebra

using BenchmarkTools

using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

function benchmark_system_my(system)
    system = GroebnerBases.change_ordering(system, :degrevlex)
    GroebnerBases.groebner([system[1]])
    @btime gb = GroebnerBases.groebner($system)
end

function benchmark_system_singular(system)
    R = parent(system[1])
    modulo = characteristic(R)
    n = nvars(R)
    ground_s = Singular.N_ZpField(modulo)
    R_s, _ = Singular.PolynomialRing(ground_s, ["x$i" for i in 1:n], ordering=:degrevlex)

    system_s = map(
        f -> AbstractAlgebra.change_base_ring(
                    ground_s,
                    AbstractAlgebra.map_coefficients(c -> ground_s(c.d), f),
                    parent=R_s),
        system)

    ideal_s = Singular.Ideal(R_s, system_s)

    Singular.std(Singular.Ideal(R_s, [system_s[1]]))

    @btime gb = Singular.std($ideal_s)
end

function run_f4_ff_degrevlex_benchmarks()
    ground = AbstractAlgebra.GF(2^31 - 1)
    systems = [
        ("cyclic 12", GroebnerBases.rootn(12, ground=ground)),
        ("cyclic 13", GroebnerBases.rootn(13, ground=ground)),
        ("katsura 9", GroebnerBases.katsura9(ground=ground)),
        ("katsura 10",GroebnerBases.katsura10(ground=ground)),
        ("noon 6"    ,GroebnerBases.noonn(6, ground=ground)),
        ("noon 7"    ,GroebnerBases.noonn(7, ground=ground))
    ]

    for (name, system) in systems
        println("$name")
        println("my")
        benchmark_system_my(system)
        println("singular")
        benchmark_system_singular(system)
    end
end
