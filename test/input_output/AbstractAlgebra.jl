using AbstractAlgebra, Test
import Primes

@testset "AbstractAlgebra.jl, univariate" begin
    R, x = polynomial_ring(GF(2^62 + 135), "x")
    @test Groebner.groebner([R(2)]) == [R(1)]
    @test Groebner.groebner([R(0), R(0)]) == [R(0)]
    @test Groebner.groebner([R(0), R(3), R(0)]) == [R(1)]

    for ground in [GF(2^31 - 1), GF(2^62 + 135), QQ]
        for gb_ord in [Groebner.Lex(), Groebner.DegLex(), Groebner.DegRevLex()]
            R, x = polynomial_ring(ground, "x")
            @test Groebner.groebner([x^2 - 4, x + 2], ordering=gb_ord) == [x + 2]
        end
    end
end

@testset "AbstractAlgebra.jl, input-output" begin
    # R, (x, y) = AbstractAlgebra.GF(Primes.nextprime(BigInt(2)^100))["x", "y"]
    # @test_throws DomainError Groebner.groebner([x, y])

    R, (x, y) = AbstractAlgebra.ZZ["x", "y"]
    @test_throws DomainError Groebner.groebner([x, y])

    R_, (a) = AbstractAlgebra.QQ["a"]
    R, (x, y) = R_["x", "y"]
    @test_throws DomainError Groebner.groebner([x, y])

    # R_, (a) = AbstractAlgebra.QQ["a"]
    # R, (x, y) = AbstractAlgebra.FractionField(R_)["x", "y"]
    # @test_throws DomainError Groebner.groebner([x, y])

    aa_orderings_to_test = [:lex, :degrevlex, :deglex]
    aa_grounds_to_test = [
        AbstractAlgebra.GF(2^62 + 135),
        AbstractAlgebra.GF(2^31 - 1),
        AbstractAlgebra.GF(17),
        AbstractAlgebra.QQ
    ]

    for aaord in aa_orderings_to_test
        R, (x, y) = AbstractAlgebra.polynomial_ring(
            AbstractAlgebra.GF(5),
            ["x", "y"],
            internal_ordering=aaord
        )
        fs = [x + 3y, R(0)]
        for gbord in [
            Groebner.Lex(x, y),
            Groebner.DegLex(x, y),
            Groebner.DegRevLex(x, y),
            Groebner.DegLex(x) * Groebner.DegRevLex(y)
        ]
            gb = Groebner.groebner(fs)
            @test gb == [x + 3y]
            @test parent(gb[1]) == R
        end
    end

    for ord in aa_orderings_to_test
        for ground in aa_grounds_to_test
            R, (x,) = polynomial_ring(ground, ["x"], internal_ordering=ord)
            @test parent(first(Groebner.groebner([x]))) == R

            R, (x, y) = polynomial_ring(ground, ["x", "y"], internal_ordering=ord)
            fs = [x^2 * y + 3, (2^31 - 5) * x - (2^31 - 4) * y]
            gb = Groebner.groebner(fs)
            @test parent(gb[1]) == R
            @test all(
                in([
                    y^3 + ground(13835057990857654347) // ground(4611686001247518736),
                    x - ground(2147483644) // ground(2147483643) * y
                ]),
                gb
            )

            @test Groebner.isgroebner(gb)
            @test all(iszero, Groebner.normalform(gb, gb))
        end
    end

    # Test for different Groebner.jl orderings
    aa_ord = :lex
    for ground in aa_grounds_to_test
        R, (x,) = polynomial_ring(ground, ["x"], internal_ordering=aa_ord)
        for case in [
            Groebner.Lex(),
            Groebner.DegLex(),
            Groebner.DegRevLex(),
            Groebner.Lex(x),
            Groebner.DegLex(x)
        ]
            gb = Groebner.groebner([x^2], ordering=case)
            @test parent(first(gb)) == R
            @test gb == [x^2]
        end

        R, (x, y) = polynomial_ring(ground, ["x", "y"], internal_ordering=aa_ord)
        fs = [x^2 + 3, y - 1]
        for case in [
            Groebner.Lex(),
            Groebner.DegLex(),
            Groebner.DegRevLex(),
            Groebner.Lex(x, y),
            Groebner.Lex(y, x)
        ]
            gb = Groebner.groebner(fs, ordering=case)
            @test parent(first(gb)) == R
            @test all(in(fs), gb)
        end
    end
end
