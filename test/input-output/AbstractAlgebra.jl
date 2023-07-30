using AbstractAlgebra
import Primes

@testset "AbstractAlgebra.jl, univariate" begin
    R, x = PolynomialRing(GF(2^62 + 135), "x")
    @test Groebner.groebner([R(2)]) == [R(1)]
    @test Groebner.groebner([R(0), R(0)]) == [R(0)]
    @test Groebner.groebner([R(0), R(3), R(0)]) == [R(1)]

    for ground in [GF(2^31 - 1), GF(2^62 + 135), QQ]
        for gb_ord in [Groebner.Lex(), Groebner.DegLex(), Groebner.DegRevLex()]
            R, x = PolynomialRing(ground, "x")
            @test Groebner.groebner([x^2 - 4, x + 2], ordering=gb_ord) == [x + 2]
        end
    end
end

@testset "AbstractAlgebra.jl, input-output" begin
    R, (x, y) = AbstractAlgebra.GF(Primes.nextprime(BigInt(2)^100))["x", "y"]
    @test_throws DomainError Groebner.groebner([x, y])

    aa_orderings_to_test = [:lex, :degrevlex, :deglex]
    aa_grounds_to_test = [
        AbstractAlgebra.GF(2^62 + 135),
        AbstractAlgebra.GF(2^31 - 1),
        AbstractAlgebra.GF(17),
        AbstractAlgebra.QQ
    ]

    for ord in aa_orderings_to_test
        for ground in aa_grounds_to_test
            R, (x,) = PolynomialRing(ground, ["x"], ordering=ord)
            @test parent(first(Groebner.groebner([x]))) == R

            R, (x, y) = PolynomialRing(ground, ["x", "y"], ordering=ord)
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
        R, (x,) = PolynomialRing(ground, ["x"], ordering=aa_ord)
        for case in [
            (gb_ord=Groebner.Lex(), same_parent=true),
            (gb_ord=Groebner.DegLex(), same_parent=false),
            (gb_ord=Groebner.DegRevLex(), same_parent=false),
            (gb_ord=Groebner.Lex(x), same_parent=true),
            (gb_ord=Groebner.DegLex(x), same_parent=false)
        ]
            gb_ord = case.gb_ord
            same_parent = case.same_parent
            gb = Groebner.groebner([x^2], ordering=gb_ord)
            if same_parent
                @test parent(first(gb)) == R
                @test gb == [x^2]
            else
                @test repr(gb[1]) == "x^2"
            end
        end

        R, (x, y) = PolynomialRing(ground, ["x", "y"], ordering=aa_ord)
        fs = [x^2 + 3, y - 1]
        for case in [
            (gb_ord=Groebner.Lex(), same_parent=true),
            (gb_ord=Groebner.DegLex(), same_parent=false),
            (gb_ord=Groebner.DegRevLex(), same_parent=false),
            (gb_ord=Groebner.Lex(x, y), same_parent=true),
            (gb_ord=Groebner.Lex(y, x), same_parent=true)
        ]
            gb_ord = case.gb_ord
            same_parent = case.same_parent
            gb = Groebner.groebner(fs, ordering=gb_ord)
            if same_parent
                @test parent(first(gb)) == R
                @test all(in(fs), gb)
            else
                @test repr(x^2 + 3) in map(repr, gb)
                @test repr(y - 1) in map(repr, gb)
            end
        end
    end
end
