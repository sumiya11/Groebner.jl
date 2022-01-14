
import Singular
import AbstractAlgebra

using FastGroebner

using BenchmarkTools

using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

function benchmark_system_my(system)
    system = FastGroebner.change_ordering(system, :degrevlex)
    FastGroebner.groebner([system[1]])
    @btime gb = FastGroebner.groebner($system, reduced=false)
    # println("length = $(length(FastGroebner.groebner(system)))")
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

    Singular.std(Singular.Ideal(R_s, [system_s[1]]))

    @btime gb = Singular.std($ideal_s, complete_reduction=false)
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 12", FastGroebner.rootn(12, ground=ground)),
        ("cyclic 13", FastGroebner.rootn(13, ground=ground)),
        ("katsura 9", FastGroebner.katsura9(ground=ground)),
        ("katsura 10",FastGroebner.katsura10(ground=ground)),
        ("noon 6"    ,FastGroebner.noonn(6, ground=ground)),
        ("noon 7"    ,FastGroebner.noonn(7, ground=ground))
    ]

    for (name, system) in systems
        println("$name")
        println("my")
        benchmark_system_my(system)
        println("singular")
        benchmark_system_singular(system)
    end
end

ground = AbstractAlgebra.GF(2^31 - 1)
run_f4_ff_degrevlex_benchmarks(ground)

#=

     readme_f4_vs_Singular.jl                                                                         cyclic 12                                                                                         my                                                                                                  834.491 ms (342159 allocations: 71.10 MiB)                                                      singular                                                                                            32.055 ms (2 allocations: 48 bytes)                                                             cyclic 13                                                                                         my                                                                                                  3.517 s (868686 allocations: 188.93 MiB)                                                        singular                                                                                            103.037 ms (2 allocations: 48 bytes)                                                            katsura 9                                                                                         my                                                                                                  717.322 ms (83858 allocations: 23.48 MiB)                                                       singular                                                                                            2.751 s (2 allocations: 48 bytes)                                                               katsura 10                                                                                        my                                                                                                  5.621 s (252687 allocations: 78.50 MiB)                                                         singular                                                                                            25.076 s (2 allocations: 48 bytes)                                                              noon 6                                                                                            my                                                                                                  91.861 ms (103667 allocations: 22.10 MiB)                                                       singular                                                                                            102.661 ms (2 allocations: 48 bytes)                                                            noon 7                                                                                            my                                                                                                  691.087 ms (501602 allocations: 104.85 MiB)                                                     singular                                                                                            896.953 ms (2 allocations: 48 bytes)
     
=#
