using DynamicPolynomials, Test

@testset "DynamicPolynomials.jl input-output" begin
    @polyvar x y

    fs = [zero(x)]
    @test Groebner.groebner(fs) == [zero(x)]

    fs = [zero(x), zero(x), zero(x)]
    @test Groebner.groebner(fs) == [zero(x)]

    fs = [zero(x), 3 * one(x), zero(x)]
    @test Groebner.groebner(fs) == [one(x)]

    fs = [2x + 3, 4y + 5]
    @test Groebner.groebner(fs) == [y + 5 // 4, x + 3 // 2]

    fs = [UInt16(2)x + UInt16(3), UInt16(2)y]
    @test Groebner.groebner(fs) == [y, x + 3 // 2]

    fs = [BigInt(1)x + BigInt(20), y]
    @test Groebner.groebner(fs) == [y, x + 20]

    fs = [x + BigInt(20) // BigInt(3), 4y + 5]
    @test Groebner.groebner(fs) == [y + 5 // 4, x + BigInt(20) // BigInt(3)]

    fs = [BigInt(34343343433)x * y^2 + 3431234567833, 3434343434x * y - 342343242342]
    @test Groebner.isgroebner(Groebner.groebner(fs))
    @test all(iszero, Groebner.normalform(Groebner.groebner(fs), fs))

    @polyvar x1 x2 x3
    noon3 = [
        10 * x1 * x2^2 + 10 * x1 * x3^2 - 11 * x1 + 10
        10 * x1^2 * x2 + 10 * x2 * x3^2 - 11 * x2 + 10
        10 * x1^2 * x3 + 10 * x2^2 * x3 - 11 * x3 + 10
    ]
    @test Groebner.groebner(noon3) == [
        x1 * x2^2 + x1 * x3^2 - 11 // 10 * x1 + 1,
        x1^2 * x3 + x2^2 * x3 - 11 // 10 * x3 + 1,
        x1^2 * x2 + x2 * x3^2 - 11 // 10 * x2 + 1,
        x2^2 * x3^2 + 11 // 20 * x1^2 - 11 // 20 * x2^2 - 11 // 20 * x3^2 - 1 // 2 * x1 +
        1 // 2 * x2 +
        1 // 2 * x3,
        x2^3 * x3 - x2 * x3^3 + x2 - x3,
        x2^4 - x3^4 - 20 // 11 * x1 * x3^2 - 10 // 11 * x2^3 + 10 // 11 * x2^2 * x3 -
        10 // 11 * x2 * x3^2 + 10 // 11 * x3^3 - 11 // 10 * x2^2 +
        11 // 10 * x3^2 +
        x1 +
        x2 - x3 - 10 // 11,
        x1^4 - x3^4 - 10 // 11 * x1^3 - 10 // 11 * x1 * x3^2 - 10 // 11 * x2^2 * x3 -
        20 // 11 * x2 * x3^2 + 10 // 11 * x3^3 - 11 // 10 * x1^2 +
        11 // 10 * x3^2 +
        x1 +
        x2 - 20 // 11,
        x3^5 + 20 // 11 * x1 * x3^3 + 20 // 11 * x2 * x3^3 - 10 // 11 * x3^4 - 33 // 20 * x3^3 +
        1 // 2 * x1^2 - 3 // 2 * x1 * x3 + 1 // 2 * x2^2 - 3 // 2 * x2 * x3 + x3^2 - 5 // 11 * x1 -
        5 // 11 * x2 + 6331 // 2200 * x3 - 11 // 20,
        x2 * x3^4 - 11 // 20 * x2^3 - 11 // 10 * x2 * x3^2 - 1 // 2 * x1 * x2 + 1 // 2 * x2^2 -
        1 // 2 * x2 * x3 +
        x3^2 +
        121 // 200 * x2 - 11 // 20,
        x1 * x3^4 - 11 // 20 * x1^3 - 11 // 10 * x1 * x3^2 + 1 // 2 * x1^2 - 1 // 2 * x1 * x2 -
        1 // 2 * x1 * x3 +
        x3^2 +
        121 // 200 * x1 - 11 // 20,
        x1 * x2 * x3^3 - 11 // 20 * x1 * x2 * x3 - 1 // 2 * x1 * x2 +
        1 // 2 * x1 * x3 +
        1 // 2 * x2 * x3
    ]

    # Test for different Groebner.jl orderings
    @polyvar x
    for gb_ord in [
        Groebner.Lex(),
        Groebner.DegLex(),
        Groebner.DegRevLex(),
        Groebner.Lex(x),
        Groebner.DegRevLex(x)
    ]
        gb = Groebner.groebner([x^2 + 1], ordering=gb_ord)
        @test gb == [x^2 + 1]
    end

    @polyvar x y
    fs = [x^2 + 3, y - 1]
    gb = Groebner.groebner(fs)
    for case in [
        (gb_ord=Groebner.Lex(), result=[y - 1, x^2 + 3]),
        (gb_ord=Groebner.DegLex(y, x), result=[y - 1, x^2 + 3]),
        (gb_ord=Groebner.DegRevLex(y, x), result=[y - 1, x^2 + 3]),
        (gb_ord=Groebner.Lex(x, y), result=[y - 1, x^2 + 3]),
        (gb_ord=Groebner.Lex(y, x), result=[x^2 + 3, y - 1])
    ]
        gb_ord = case.gb_ord
        result = case.result
        gb = Groebner.groebner(fs, ordering=gb_ord)
        @test gb == result
    end

    @polyvar x y z
    Groebner.groebner([x, y, z^2], ordering=Groebner.Lex(x, y, z)) == [z^2, y, x]
    Groebner.groebner([x, y, z^2], ordering=Groebner.Lex(y, x, z)) == [z^2, x, y]
    Groebner.groebner([x, y, z^2], ordering=Groebner.Lex(z, x, y)) == [y, x, z^2]

    # Test for different DynamicPolynomials.jl orderings
    @polyvar x y monomial_order = LexOrder
    F = [2x^2 * y + 3x, 4x * y^2 + 5y^3]
    gb = Groebner.groebner(F, ordering=Groebner.DegLex())
    @test gb == [(5 // 4) * x * y + x^2, y^3 + (24 // 25) * x, (-6 // 5) * x + x * y^2]
    gb = Groebner.groebner(F, ordering=Groebner.DegRevLex())
    @test gb == [(5 // 4) * x * y + x^2, y^3 + (24 // 25) * x, (-6 // 5) * x + x * y^2]
    gb = Groebner.groebner(F, ordering=Groebner.Lex())
    @test gb == [(-6 // 5) * y^3 + y^5, (25 // 24) * y^3 + x]
    gb = Groebner.groebner(F, ordering=Groebner.Lex(y, x))
    @test gb == [(-15 // 8) * x + x^3, x * y + (4 // 5) * x^2, (24 // 25) * x + y^3]

    # It is possible to over-specify the ordering
    @polyvar x y z
    @test groebner([x, y], ordering=Lex(x, y, z)) == [y, x]
    @test groebner([x, y], ordering=Lex(y, x, z)) == [x, y]

    @polyvar x y monomial_order = Graded{Reverse{InverseLexOrder}}
    F = [2x^2 * y + 3x, 4x * y^2 + 5y^3]
    gb = Groebner.groebner(F)
    @test gb == [(5 // 4) * x * y + x^2, y^3 + (24 // 25) * x, (-6 // 5) * x + x * y^2]
    gb = Groebner.groebner(F, ordering=Groebner.Lex())
    @test gb == [(-6 // 5) * y^3 + y^5, (25 // 24) * y^3 + x]

    for dp_ord in [LexOrder, Graded{LexOrder}, Graded{Reverse{InverseLexOrder}}]
        @polyvar x y monomial_order = dp_ord
        F = [2x^2 * y + 3x, 4x * y^2 + 5y^3]
        for gb_ord in [
            Groebner.Lex(),
            Groebner.DegLex(),
            Groebner.DegRevLex(),
            Groebner.Lex(x, y),
            Groebner.Lex(y, x)
        ]
            gb = Groebner.groebner(F, ordering=gb_ord)
            @test Groebner.isgroebner(gb, ordering=gb_ord)
        end
    end

    @polyvar x y monomial_order = InverseLexOrder
    F = [2x^2 * y + 3x, 4x * y^2 + 5y^3]
    @test_throws DomainError Groebner.groebner(F)
end
