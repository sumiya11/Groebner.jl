using DynamicPolynomials

@testset "DynamicPolynomials.jl input-output" begin
    @polyvar x y

    fs = [2x + 3, 4y + 5]
    @test Groebner.groebner(fs) == [4y + 5, 2x + 3]

    fs = [UInt16(2)x + UInt16(3), UInt16(2)y]
    @test Groebner.groebner(fs) == [y, UInt16(2)x + UInt16(3)]

    fs = [BigInt(1)x + BigInt(20), y]
    @test Groebner.groebner(fs) == [y, BigInt(1)x + BigInt(20)]

    fs = [x + BigInt(20) // BigInt(3), 4y + 5]
    @test Groebner.groebner(fs) == [y + 5 // 4, x + BigInt(20) // BigInt(3)]

    fs = [34343343433x * y^2 + 3431234567833, 3434343434x * y - 342343242342]
    @test_throws DomainError Groebner.groebner(fs)

    fs = [BigInt(34343343433)x * y^2 + 3431234567833, 3434343434x * y - 342343242342]
    @test Groebner.isgroebner(Groebner.groebner(fs))
    @test all(iszero, Groebner.normalform(Groebner.groebner(fs), fs))

    @polyvar x1 x2 x3
    noon3 = [
        10 * x1 * x2^2 + 10 * x1 * x3^2 - 11 * x1 + 10 // 1
        10 * x1^2 * x2 + 10 * x2 * x3^2 - 11 * x2 + 10 // 1
        10 * x1^2 * x3 + 10 * x2^2 * x3 - 11 * x3 + 10 // 1
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
        x3^5 + 20 // 11 * x1 * x3^3 + 20 // 11 * x2 * x3^3 - 10 // 11 * x3^4 -
        33 // 20 * x3^3 + 1 // 2 * x1^2 - 3 // 2 * x1 * x3 + 1 // 2 * x2^2 -
        3 // 2 * x2 * x3 + x3^2 - 5 // 11 * x1 - 5 // 11 * x2 + 6331 // 2200 * x3 -
        11 // 20,
        x2 * x3^4 - 11 // 20 * x2^3 - 11 // 10 * x2 * x3^2 - 1 // 2 * x1 * x2 +
        1 // 2 * x2^2 - 1 // 2 * x2 * x3 +
        x3^2 +
        121 // 200 * x2 - 11 // 20,
        x1 * x3^4 - 11 // 20 * x1^3 - 11 // 10 * x1 * x3^2 + 1 // 2 * x1^2 -
        1 // 2 * x1 * x2 - 1 // 2 * x1 * x3 +
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
end
