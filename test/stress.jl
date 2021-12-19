
using .GroebnerBases: f4, is_groebner

function generate_set(nvariables, exps, nterms, npolys;
                            ground=GF(2^31 - 1), ordering=:lex)

    R, _ = PolynomialRing(
        ground,
        ["x$i" for i in 1:nvariables],
        ordering=ordering
    )

    return [
        rand(R, exps, nterms)
        for _ in 1:rand(npolys)
    ]
end


@testset "F4 random stress tests" begin


    nvariables = [2, 3]
    exps       = [1:2, 2:4, 1:4]
    nterms     = [1:1, 1:2, 2:4]
    npolys     = [1:1, 1:3, 2:4]
    grounds    = [GF(1031), GF(2^31 - 1)]
    orderings  = [:degrevlex]

    p = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings)))

    @debug "producing $p small tests for f4"

    for n in nvariables
        for e in exps
            for nt in nterms
                for np in npolys
                    for gr in grounds
                        for ord in orderings
                            set = generate_set(
                                n, e, nt, np, ground=gr, ordering=ord
                            )

                            gb = f4(set)
                            @test is_groebner(gb, initial_gens=set)

                            if !is_groebner(gb, initial_gens=set)
                                @error "BEDA" n e nt np gr ord
                            end
                        end
                    end
                end
            end
        end
    end


end
