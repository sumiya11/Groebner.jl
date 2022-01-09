
import Random

function generate_set(nvariables, exps, nterms, npolys, csz, rng;
                            ground=GF(2^31 - 1), ordering=:lex)

    R, _ = PolynomialRing(
        ground,
        ["x$i" for i in 1:nvariables],
        ordering=ordering
    )

    return filter!(!iszero, [
        map_coefficients(
            c -> ground(data(c) % csz),
            rand(rng, R, exps, nterms)
        )
        for _ in 1:rand(rng, npolys)
    ])
end


@testset "f4 random stress tests" begin

    rng = Random.MersenneTwister(5)

    nvariables = [2, 3, 4, 5]
    exps       = [1:2, 2:4, 1:4, 5:6]
    nterms     = [1:1, 1:2, 2:4, 5:6]
    npolys     = [1:1, 1:3, 2:4, 20:21]
    grounds    = [GF(1031), GF(2^31 - 1)]
    coeffssize = [3]
    orderings  = [:degrevlex]

    p = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings, coeffssize)))

    @info "producing $p small tests for f4"

    for n in nvariables
        for e in exps
            for nt in nterms
                for np in npolys
                    for gr in grounds
                        for ord in orderings
                            for csz in coeffssize
                                set = generate_set(
                                    n, e, nt, np, csz, rng, ground=gr, ordering=ord
                                )
                                isempty(set) && continue

                                gb = GroebnerBases.groebner(set)
                                @test GroebnerBases.isgroebner(gb, initial_gens=set)

                                if !GroebnerBases.isgroebner(gb, initial_gens=set)
                                    @error "Beda!" n e nt np gr ord
                                    println("Rng ", rng)
                                    println("Set ", set)
                                    println("Gb ", gb)
                                    error("exit tests")
                                end
                            end
                        end
                    end
                end
            end
        end
    end


end
