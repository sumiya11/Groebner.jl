import Nemo
import Primes

@testset "Nemo.jl, univariate" begin
    for ff in [Nemo.Native.GF(2^62 + 135), Nemo.GF(2^62 + 135)]
        R, x = polynomial_ring(ff, "x")
        @test Groebner.groebner([R(2)]) == [R(1)]
        @test Groebner.groebner([R(0), R(0)]) == [R(0)]
        @test Groebner.groebner([R(0), R(3), R(0)]) == [R(1)]
        @test Groebner.groebner([R(0), R(3), R(0)]) == [R(1)]
        @test Groebner.groebner([x + 100]) == [x + 100]
    end

    for ff in [
        Nemo.Native.GF(2^31 - 1),
        Nemo.Native.GF(2^62 + 135),
        Nemo.GF(2^62 + 135),
        Nemo.GF(Nemo.ZZRingElem(2^62 + 135)),
        Nemo.QQ
    ]
        for gb_ord in [Groebner.Lex(), Groebner.DegLex(), Groebner.DegRevLex()]
            R, x = polynomial_ring(ff, "x")
            @test Groebner.groebner([x^2 - 4, x + 2], ordering=gb_ord) == [x + 2]
        end
    end
end

@testset "Nemo.jl, input-output" begin
    R, (x, y) = Nemo.GF(Nemo.ZZRingElem(Primes.nextprime(BigInt(2)^100)))["x", "y"]
    @test_throws DomainError Groebner.groebner([x, y])

    R, (x, y) = Nemo.GF(2, 2)["x", "y"]
    @test_throws DomainError Groebner.groebner([x, y])

    R, (x, y) = Nemo.ZZ["x", "y"]
    @test_throws DomainError Groebner.groebner([x, y])

    R_, (a) = Nemo.QQ["a"]
    R, (x, y) = R_["x", "y"]
    @test_throws DomainError Groebner.groebner([x, y])

    # R_, (a) = Nemo.QQ["a"]
    # R, (x, y) = AbstractAlgebra.FractionField(R_)["x", "y"]
    # @test_throws DomainError Groebner.groebner([x, y])

    nemo_orderings_to_test = [:lex, :deglex, :degrevlex]
    nemo_grounds_to_test = [
        Nemo.Native.GF(2^62 + 135),
        Nemo.Native.GF(2^31 - 1),
        Nemo.Native.GF(17),
        Nemo.GF(2^62 + 135),
        Nemo.GF(2^31 - 1),
        Nemo.QQ
    ]

    for ord in nemo_orderings_to_test
        for ground in nemo_grounds_to_test
            R, x = polynomial_ring(ground, "x")
            gb = Groebner.groebner([(x - 1) * (x + 8), (x + 8) * (x + 10)])
            @test gb == [(x + 8)]

            R, (x, y) = polynomial_ring(ground, ["x", "y"], ordering=ord)
            fs = [x^2 * y + 3, (2^31 - 5) * x - (2^31 - 4) * y]
            gb = Groebner.groebner(fs)
            @test parent(gb[1]) == R
            @test Groebner.isgroebner(gb)
        end
    end

    # Test for different Groebner.jl orderings
    for nemo_ord in [:lex, :deglex, :degrevlex]
        for ground in nemo_grounds_to_test
            R, (x,) = polynomial_ring(ground, ["x"], ordering=nemo_ord)
            for gb_ord in [
                Groebner.Lex(),
                Groebner.DegLex(),
                Groebner.DegRevLex(),
                Groebner.Lex(x),
                Groebner.DegLex(x)
            ]
                gb = Groebner.groebner([x^2], ordering=gb_ord)
                @test parent(first(gb)) == R
                @test gb == [x^2]
            end

            R, (x, y) = polynomial_ring(ground, ["x", "y"], ordering=nemo_ord)
            fs = [x^2 + 3, y - 1]
            for gb_ord in [
                Groebner.Lex(),
                Groebner.DegLex(),
                Groebner.DegRevLex(),
                Groebner.Lex(x, y),
                Groebner.Lex(y, x)
            ]
                gb = Groebner.groebner(fs, ordering=gb_ord)
                @test parent(first(gb)) == R
                @test all(in(fs), gb) && all(in(gb), fs)
            end
        end
    end
end
