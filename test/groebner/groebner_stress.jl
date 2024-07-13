import Random
using AbstractAlgebra

function test_params(
    rng,
    nvariables,
    maxdegs,
    nterms,
    npolys,
    grounds,
    coeffssize,
    orderings,
    linalgs,
    monoms,
    homogenizes
)
    boot = 1
    for _ in 1:boot
        for nv in nvariables
            for md in maxdegs
                for nt in nterms
                    for np in npolys
                        for k in grounds
                            for ord in orderings
                                for cf in coeffssize
                                    set = Groebner.Examples.random_generating_set(
                                        rng,
                                        k,
                                        ord,
                                        nv,
                                        md,
                                        nt,
                                        np,
                                        cf
                                    )
                                    isempty(set) && continue

                                    for linalg in linalgs
                                        for monom in monoms
                                            for homogenize in homogenizes
                                                try
                                                    gb = Groebner.groebner(
                                                        set,
                                                        linalg=linalg,
                                                        monoms=monom,
                                                        homogenize=homogenize
                                                    )
                                                    flag = Groebner.isgroebner(gb)
                                                    if !flag
                                                        @error "Beda!" nv md nt np k ord monom
                                                        println("Rng:\n", rng)
                                                        println("Set:\n", set)
                                                        println("Gb:\n", gb)
                                                    end
                                                    @test flag
                                                catch err
                                                    @error "Beda!" nv md nt np k ord monom
                                                    println(err)
                                                    println("Rng:\n", rng)
                                                    println("Set:\n", set)
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
        end
    end
end

@testset "groebner random stress tests" begin
    rng = Random.MersenneTwister(42)

    nvariables = [2, 3]
    maxdegs    = [2, 4]
    nterms     = [1, 2, 4]
    npolys     = [1, 4, 100]
    grounds    = [GF(1031), GF(2^50 + 55), AbstractAlgebra.QQ]
    coeffssize = [3, 1000, 2^31 - 1, BigInt(2)^100]
    orderings  = [:degrevlex, :lex, :deglex]
    linalgs    = [:deterministic, :randomized]
    monoms     = [:auto, :dense, :packed]
    homogenize = [:yes, :auto]
    p          = prod(map(length, (nvariables, maxdegs, nterms, npolys, grounds, orderings, coeffssize, linalgs, monoms, homogenize)))
    @info "Producing $p small random tests for groebner. This may take a minute"

    test_params(
        rng,
        nvariables,
        maxdegs,
        nterms,
        npolys,
        grounds,
        coeffssize,
        orderings,
        linalgs,
        monoms,
        homogenize
    )
end
