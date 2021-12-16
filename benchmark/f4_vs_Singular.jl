
import GroebnerBases: katsura6, change_ordering
import AbstractAlgebra
import Singular

using BenchmarkTools


function run_ours(system)
    @benchmark gb = GroebnerBases.f4(system)
end

function run_singular(system)
    @benchmark gb = Singular.std(system)
end

function benchmark_katsura()
    modulo = 2^31 - 1
    system = change_ordering(katsura6(ground=AbstractAlgebra.GF(2^31-1), :degrevlex))

    ground_s = Singular.N_ZpField(2^31 - 1)
    R_s, _ = Singular.PolynomialRing(ground_s, ["x$i" for i in 1:7], ordering=:degrevlex)

    system_s = map(
        f -> change_base_ring(ground_s, f, parent=R_s),
        system)

    ideal_s = Singular.Ideal(R_s, system_s)

    run_ours(system)
    run_singular(ideal_s)
end
