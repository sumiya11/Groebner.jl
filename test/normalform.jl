using AbstractAlgebra, Random

@testset "normalform" begin
    R, x = polynomial_ring(GF(2^31 - 1), "x")

    @test Groebner.normalform([x], R(0)) == R(0)
    @test Groebner.normalform([x], R(1)) == R(1)
    @test Groebner.normalform([x], [R(0)]) == [R(0)]
    @test Groebner.normalform([x], [R(1)]) == [R(1)]
    @test Groebner.normalform([x], [x, R(5), x + 1]) == [R(0), R(5), R(1)]

    R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"])

    @test Groebner.normalform([x], R(0)) == R(0)
    @test Groebner.normalform([x], R(1)) == R(1)
    @test Groebner.normalform([x], [R(0)]) == [R(0)]
    @test Groebner.normalform([x], [R(1)]) == [R(1)]

    Gs = [[x], [x, y], [x, y, z]]

    @test Groebner.normalform(Gs[1], x) == 0
    @test Groebner.normalform(Gs[1], y) == y
    @test Groebner.normalform(Gs[3], x + y^2) == 0
    @test Groebner.normalform(Gs[2], 5x^2 - y^3 - 1) == -1
    @test Groebner.normalform(Gs[2], x + y + 3z) == 3z

    G = [x^2 + y * x + 1, y^2 - y]

    @test Groebner.normalform(G, x + y + z) == x + y + z
    @test Groebner.normalform(G, x^2 + y^2) == y - y * x - 1
    @test Groebner.normalform(G, y^3) == Groebner.normalform(G, y^4) == y

    R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"], internal_ordering=:degrevlex)
    G = [x^2 + y, y^2 + x]
    @test Groebner.normalform(G, x^2 + y^2) == -x - y

    @test_nowarn Groebner.normalform([x, x * y], x)
end

@testset "normalform many variables" begin
    R, xs = polynomial_ring(QQ, ["x$i" for i in 1:63])
    GB = Groebner.groebner(xs)
    @test GB == reverse(xs)
    @test Groebner.normalform(GB, xs[1]) == R(0)

    R, xs = polynomial_ring(QQ, ["x$i" for i in 1:220])
    GB = Groebner.groebner(xs)
    @test GB == reverse(xs)
    @test Groebner.normalform(GB, xs[1]) == R(0)

    R, xs = polynomial_ring(QQ, ["x$i" for i in 1:1025])
    GB = Groebner.groebner(xs)
    @test GB == reverse(xs)
    @test Groebner.normalform(GB, xs) == zeros(R, length(xs))
end

@testset "normalform orderings" begin
    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])

    @test Groebner.normalform([x + y], [x, y]) == [-y, y]
    @test Groebner.normalform([x + y^2], [x, y], ordering=Groebner.DegLex()) == [x, y]
    @test Groebner.normalform([x + y^2], [x, y], ordering=Groebner.DegRevLex()) == [x, y]
    @test Groebner.normalform([x + y], [x, y], ordering=Groebner.Lex(y, x, z)) == [x, -x]
end

@testset "normalform of an array" begin
    for field in [GF(17), GF(2^31 - 1), QQ]
        R, (x, y) = polynomial_ring(field, ["x", "y"])
        gb = [x, y]
        @test Groebner.normalform(gb, [x, y + 1]) == [R(0), R(1)]
        @test Groebner.normalform(gb, [y + 1, x]) == [R(1), R(0)]
        @test Groebner.normalform(gb, [x, 3y, 4x + 5y, 6x * y, 7x + 1, y + 2, x^2 + y^2 + 3]) ==
              [R(0), R(0), R(0), R(0), R(1), R(2), R(3)]

        gb = [x]
        @test Groebner.normalform(gb, [x, 3y, 4x + 5y, 6x * y, 7x + 1, y + 2, x^2 + y^2 + 3]) ==
              [R(0), 3y, 5y, R(0), R(1), y + 2, y^2 + 3]

        gb = [x^2 + y^2, y^3 - y]
        @test Groebner.normalform(
            gb,
            [R(1), R(5), R(16), 3x, y, 4x * y + 2, y^2, x * y^2 - 8x + y]
        ) == [R(1), R(5), R(16), 3x, y, 4x * y + 2, y^2, x * y^2 - 8x + y]

        gb = [x, R(0), y, R(0)]
        @test Groebner.normalform(gb, [R(0)]) == [R(0)]
        @test Groebner.normalform(gb, [R(0), R(1), R(0)]) == [R(0), R(1), R(0)]
        @test Groebner.normalform([R(0), R(0)], [R(0), R(1), R(0)]) == [R(0), R(1), R(0)]
    end
end

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
                                            Groebner.normalform(gb, f, ordering=ord_groebner)
                                        )
                                        @test all(
                                            iszero,
                                            Groebner.normalform(gb, gb, ordering=ord_groebner)
                                        )

                                        isempty(set2) && continue
                                        f = rand(set2)
                                        f_nf = Groebner.normalform(gb, f, ordering=ord_groebner)
                                        @test iszero(
                                            Groebner.normalform(
                                                gb,
                                                f - f_nf,
                                                ordering=ord_groebner
                                            )
                                        )
                                    catch
                                        @error "Beda!" n e nt np gr ord ord_groebner
                                        println(err)
                                        println("Rng:\n", rng)
                                        println("Set:\n", set1)
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
            (nvariables, exps, nterms, npolys, grounds, orderings, orderings_groebner, coeffssize)
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
