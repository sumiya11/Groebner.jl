
import Random
using AbstractAlgebra

function test_params(rng, nvariables, exps, nterms,
                        npolys, grounds, coeffssize, orderings, linalgs)

    for n in nvariables
        for e in exps
            for nt in nterms
                for np in npolys
                    for gr in grounds
                        for ord in orderings
                            for csz in coeffssize
                                for linalg in linalgs
                                    set = Groebner.generate_set(
                                        n, e, nt, np, csz, rng, gr, ord
                                    )
                                    filter!(!iszero, set)

                                    isempty(set) && continue
                                    # println(ord, " ", gr)
                                    # println(set)
                                    gb = Groebner.groebner(set, linalg=linalg)
                                    @test Groebner._isgroebner_reference(gb, initial_gens=set)

                                    if !Groebner._isgroebner_reference(gb, initial_gens=set)
                                        @error "Beda!" n e nt np gr ord
                                        println("Rng ", rng)
                                        println("Set ", set)
                                        println("Gb ", gb)
                                        # error("exit tests")
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


@testset "f4 random stress tests" begin

    rng = Random.MersenneTwister(42)

    # :degrevlex finite case tests
    nvariables = [2, 3, 4]
    exps       = [1:2, 2:4, 4:5]
    nterms     = [1:1, 1:2, 4:5]
    npolys     = [1:1, 1:3, 4:5]
    grounds    = [GF(1031), GF(2^31 - 1)]
    coeffssize = [3, 1000, 2^31 - 1]
    orderings  = [:degrevlex]
    linalgs    = [:exact]
    p = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings, coeffssize, linalgs)))
    @info "producing $p $(orderings[1]) tests for f4"
    test_params(rng, nvariables, exps, nterms, npolys, grounds, coeffssize, orderings, linalgs)

    # :lex finite case tests
    nvariables = [2, 3, 4]
    exps       = [1:2, 2:4, 2:3]
    nterms     = [1:1, 1:2, 2:3]
    npolys     = [1:1, 1:3, 2:3]
    grounds    = [GF(1031), GF(2^31 - 1)]
    coeffssize = [3, 1000, 2^31 - 1]
    orderings  = [:lex]
    linalgs    = [:exact]
    p = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings, coeffssize, linalgs)))
    @info "producing $p $(orderings[1]) tests for f4"
    test_params(rng, nvariables, exps, nterms, npolys, grounds, coeffssize, orderings, linalgs)

    # :deglex finite case tests
    nvariables = [2, 3, 4]
    exps       = [1:2, 2:4, 2:3]
    nterms     = [1:1, 1:2, 2:3]
    npolys     = [1:1, 1:3, 2:3]
    grounds    = [GF(1031), GF(2^31 - 1)]
    coeffssize = [3, 1000, 2^31 - 1]
    orderings  = [:deglex]
    linalgs    = [:exact]
    p = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings, coeffssize, linalgs)))
    @info "producing $p $(orderings[1]) tests for f4"
    test_params(rng, nvariables, exps, nterms, npolys, grounds, coeffssize, orderings, linalgs)

    # :degrevlex rational case tests
    nvariables = [2, 3, 4]
    exps       = [1:2, 2:4, 2:3]
    nterms     = [1:1, 1:2, 2:3]
    npolys     = [1:1, 1:3, 2:3]
    grounds    = [QQ]
    coeffssize = [3, 1000, 2^31 - 1, 9223372036854775837]
    orderings  = [:degrevlex]
    linalgs    = [:exact]
    p = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings, coeffssize, linalgs)))
    @info "producing $p $(orderings[1]) rational tests for f4"
    test_params(rng, nvariables, exps, nterms, npolys, grounds, coeffssize, orderings, linalgs)

    # :deglex rational case tests
    nvariables = [2, 3, 4]
    exps       = [1:2, 2:4, 2:3]
    nterms     = [1:1, 1:2, 2:3]
    npolys     = [1:1, 1:3, 2:3]
    grounds    = [QQ]
    coeffssize = [3, 1000, 2^31 - 1, 9223372036854775837]
    orderings  = [:deglex]
    p = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings, coeffssize)))
    @info "producing $p $(orderings[1]) rational tests for f4"
    # test_params(rng, nvariables, exps, nterms, npolys, grounds, coeffssize, orderings)

end
