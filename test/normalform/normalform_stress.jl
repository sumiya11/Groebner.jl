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
    orderings
)
    for n in nvariables
        for e in exps
            for nt in nterms
                for np in npolys
                    for gr in grounds
                        for ord in orderings
                            for csz in coeffssize
                                set = Groebner.generate_set(n, e, nt, np, csz, rng, gr, ord)
                                isempty(set) && continue

                                try
                                    gb = Groebner.groebner(set)
                                    f = rand(set)

                                    # TODO : add more tests
                                    @test Groebner.normalform(gb, f, check=false) == 0

                                    nf = Groebner.normalform(gb, gb, check=false)
                                    @test all(iszero, nf)
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

@testset "normalform random stress tests" begin
    rng = Random.MersenneTwister(42)

    nvariables = [2, 3, 4]
    exps       = [1:2, 2:4, 2:3]
    nterms     = [1:1, 1:2, 2:3]
    npolys     = [1:1, 1:3, 2:3]
    grounds    = [GF(2^31 - 1), QQ]
    coeffssize = [3, 1000, 2^31 - 1]
    orderings  = [:deglex, :lex, :degrevlex]
    p          = prod(map(length, (nvariables, exps, nterms, npolys, grounds, orderings, coeffssize)))
    @info "Producing $p tests for normal form"
    test_params_nf(rng, nvariables, exps, nterms, npolys, grounds, coeffssize, orderings)
end
