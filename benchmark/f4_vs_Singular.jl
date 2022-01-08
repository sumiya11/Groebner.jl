

# hehe
import Singular
import AbstractAlgebra

using BenchmarkTools


function benchmark_katsura()
    modulo = 2^31 - 1
    system = change_ordering(katsura6(ground=AbstractAlgebra.GF(modulo)), :degrevlex)

    ground_s = Singular.N_ZpField(modulo)
    R_s, _ = Singular.PolynomialRing(ground_s, ["x$i" for i in 1:7], ordering=:degrevlex)

    system_s = map(
        f -> AbstractAlgebra.change_base_ring(
                    ground_s,
                    AbstractAlgebra.map_coefficients(c -> ground_s(c.d), f),
                    parent=R_s),
        system)

    ideal_s = Singular.Ideal(R_s, system_s)

    GroebnerBases.f4([system[1]])
    Singular.std(Singular.Ideal(R_s, [system_s[1]]))

    @benchmark gb1 = GroebnerBases.f4($system)
    @benchmark gb2 = Singular.std($ideal_s)

end

function benchmark_root()
    modulo = 2^31 - 1
    system = change_ordering(rootn(13, ground=AbstractAlgebra.GF(modulo)), :degrevlex)

    ground_s = Singular.N_ZpField(modulo)
    R_s, _ = Singular.PolynomialRing(ground_s, ["x$i" for i in 1:13], ordering=:degrevlex)

    system_s = map(
        f -> AbstractAlgebra.change_base_ring(
                    ground_s,
                    AbstractAlgebra.map_coefficients(c -> ground_s(c.d), f),
                    parent=R_s),
        system)

    ideal_s = Singular.Ideal(R_s, system_s)

    GroebnerBases.f4([system[1]])
    Singular.std(Singular.Ideal(R_s, [system_s[1]]))

    @benchmark gb1 = GroebnerBases.f4($system)
    @benchmark gb2 = Singular.std($ideal_s)

end

benchmark_katsura()
benchmark_root()
