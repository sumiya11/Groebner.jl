
import Random
using AbstractAlgebra

function test_params_isgb(rng, nvariables, exps, nterms,
                        npolys, grounds, coeffssize, orderings)

    without(x, k) = x[1:end .!= k]

    for n in nvariables
        for e in exps
            for nt in nterms
                for np in npolys
                    for gr in grounds
                        for ord in orderings
                            for csz in coeffssize
                                set = Groebner.generate_set(
                                    n, e, nt, np, csz, rng, gr, ord
                                )
                                isempty(set) && continue
                                gb = Groebner.groebner(set)

                                @test Groebner.isgroebner(gb)
                            end
                        end
                    end
                end
            end
        end
    end
end


@testset "isgroebner random stress tests" begin

    rng = Random.MersenneTwister(5)

    # :degrevlex finite case tests
    nvariables = [2, 3, 5]
    exps       = [1:2, 2:4, 4:5]
    nterms     = [1:1, 1:2, 4:5]
    npolys     = [1:1, 1:3, 4:5]
    grounds    = [GF(1031), GF(2^31 - 1)]
    coeffssize = [3, 1000, 2^31 - 1]
    orderings  = [:degrevlex]
    p = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings, coeffssize)))
    @info "producing $p $(orderings[1]) tests for isgroebner"
    test_params_isgb(rng, nvariables, exps, nterms, npolys, grounds, coeffssize, orderings)

    # :lex finite case tests
    nvariables = [2, 3, 4]
    exps       = [1:2, 2:4, 2:3]
    nterms     = [1:1, 1:2, 2:3]
    npolys     = [1:1, 1:3, 2:3]
    grounds    = [GF(1031), GF(2^31 - 1)]
    coeffssize = [3, 1000, 2^31 - 1]
    orderings  = [:lex]
    p = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings, coeffssize)))
    @info "producing $p $(orderings[1]) tests for isgroebner"
    test_params_isgb(rng, nvariables, exps, nterms, npolys, grounds, coeffssize, orderings)

    # :deglex finite case tests
    nvariables = [2, 3, 4]
    exps       = [1:2, 2:4, 2:3]
    nterms     = [1:1, 1:2, 2:3]
    npolys     = [1:1, 1:3, 2:3]
    grounds    = [GF(1031), GF(2^31 - 1)]
    coeffssize = [3, 1000, 2^31 - 1]
    orderings  = [:deglex]
    p = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings, coeffssize)))
    @info "producing $p $(orderings[1]) tests for isgroebner"
    test_params_isgb(rng, nvariables, exps, nterms, npolys, grounds, coeffssize, orderings)


end
