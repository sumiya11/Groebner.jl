import Random
using AbstractAlgebra

function test_params(
    rng,
    nvariables,
    exps,
    nterms,
    npolys,
    grounds,
    coeffssize,
    orderings,
    linalgs
)
    for n in nvariables
        for e in exps
            for nt in nterms
                for np in npolys
                    for gr in grounds
                        for ord in orderings
                            for csz in coeffssize
                                for linalg in linalgs
                                    set = Groebner.generate_set(
                                        n,
                                        e,
                                        nt,
                                        np,
                                        csz,
                                        rng,
                                        gr,
                                        ord
                                    )
                                    filter!(!iszero, set)
                                    isempty(set) && continue

                                    try
                                        gb = Groebner.groebner(set, linalg=linalg)
                                        @test Groebner.isgroebner(gb)
                                    catch err
                                        @error "Beda!" n e nt np gr ord
                                        println(err)
                                        println("Rng:\n", rng)
                                        println("Set:\n", set)
                                        println("Gb:\n", gb)
                                        rethrow(err)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

@testset "groebner random stress tests" begin
    rng = Random.Xoshiro(42)

    nvariables = [2, 3]
    exps       = [1:2, 2:4]
    nterms     = [1:1, 1:2, 3:4]
    npolys     = [1:1, 1:3, 3:4, 100:110]
    grounds    = [GF(1031), GF(2^31 - 1), GF(2^50 + 55), AbstractAlgebra.QQ]
    coeffssize = [3, 1000, 2^31 - 1]
    orderings  = [:degrevlex, :lex, :deglex]
    linalgs    = [:deterministic, :randomized]
    p          = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings, coeffssize, linalgs)))
    @info "Producing $p random small tests for groebner"
    test_params(
        rng,
        nvariables,
        exps,
        nterms,
        npolys,
        grounds,
        coeffssize,
        orderings,
        linalgs
    )
end
