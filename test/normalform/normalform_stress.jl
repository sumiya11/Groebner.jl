import Random
using AbstractAlgebra

function test_params_nf(
    rng,
    nvariables,
    exps,
    nterms,
    npolys,
    grounds,
    coeffssize,
    orderings,
    orderings_groebner
)
    for n in nvariables
        for e in exps
            for nt in nterms
                for np in npolys
                    for gr in grounds
                        for ord in orderings
                            for ord_groebner in orderings_groebner
                                for csz in coeffssize
                                    set1 = Groebner.Examples.random_generating_set(
                                        rng,
                                        gr,
                                        ord,
                                        n,
                                        e,
                                        nt,
                                        np,
                                        csz
                                    )
                                    set2 = Groebner.Examples.random_generating_set(
                                        rng,
                                        gr,
                                        ord,
                                        n,
                                        e,
                                        nt,
                                        np,
                                        csz
                                    )
                                    isempty(set1) && continue

                                    try
                                        gb = Groebner.groebner(set1, ordering=ord_groebner)
                                        f = rand(set1)

                                        @test iszero(
                                            Groebner.normalform(
                                                gb,
                                                f,
                                                ordering=ord_groebner
                                            )
                                        )
                                        @test all(
                                            iszero,
                                            Groebner.normalform(
                                                gb,
                                                gb,
                                                ordering=ord_groebner
                                            )
                                        )

                                        isempty(set2) && continue
                                        f = rand(set2)
                                        f_nf = Groebner.normalform(
                                            gb,
                                            f,
                                            ordering=ord_groebner
                                        )
                                        @test iszero(
                                            Groebner.normalform(
                                                gb,
                                                f - f_nf,
                                                ordering=ord_groebner
                                            )
                                        )
                                    catch
                                        @error "Beda!" n e nt np gr ord
                                        println(err)
                                        println("Rng:\n", rng)
                                        println("Set:\n", set)
                                        println("Gb:\n", gb)
                                        println("Nf:\n", nf)
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

@testset "normalform random stress tests" begin
    rng = Random.MersenneTwister(42)

    nvariables = [2, 3]
    exps = [2, 4, 3]
    nterms = [2, 3]
    npolys = [2, 3]
    grounds = [GF(2^31 - 1), QQ]
    coeffssize = [3, 1000, 2^31 - 1]
    orderings = [:deglex, :lex, :degrevlex]
    orderings_groebner =
        [Groebner.InputOrdering(), Groebner.Lex(), Groebner.DegLex(), Groebner.DegRevLex()]
    p = prod(
        map(
            length,
            (
                nvariables,
                exps,
                nterms,
                npolys,
                grounds,
                orderings,
                orderings_groebner,
                coeffssize
            )
        )
    )
    @info "Producing $p tests for normal form. This may take a minute"
    test_params_nf(
        rng,
        nvariables,
        exps,
        nterms,
        npolys,
        grounds,
        coeffssize,
        orderings,
        orderings_groebner
    )
end
