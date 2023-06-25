
@testset "ff f4 certify" begin
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"], ordering=:lex)

    @test Groebner.groebner([x, y], certify=true) == [y, x]
    @test Groebner.groebner([y, x], certify=true) == [y, x]

    fs = [x^2 + y, x * y]
    @test Groebner.groebner(fs, certify=true) == [y^2, x * y, x^2 + y]
end

@testset "qq f4 certify" begin
    root = Groebner.change_ordering(Groebner.rootn(3, ground=QQ), :degrevlex)
    x1, x2, x3 = gens(parent(first(root)))
    gb = Groebner.groebner(root, certify=true)
    @test gb == [x1 + x2 + x3, x2^2 + x2 * x3 + x3^2, x3^3 - 1]

    root = Groebner.change_ordering(Groebner.rootn(6, ground=QQ), :deglex)
    x1, x2, x3, x4, x5, x6 = gens(parent(first(root)))
    gb = Groebner.groebner(root, certify=true)
    @test gb == [
        x1 + x2 + x3 + x4 + x5 + x6,
        x2^2 +
        x2 * x3 +
        x3^2 +
        x2 * x4 +
        x3 * x4 +
        x4^2 +
        x2 * x5 +
        x3 * x5 +
        x4 * x5 +
        x5^2 +
        x2 * x6 +
        x3 * x6 +
        x4 * x6 +
        x5 * x6 +
        x6^2,
        x3^3 +
        x3^2 * x4 +
        x3 * x4^2 +
        x4^3 +
        x3^2 * x5 +
        x3 * x4 * x5 +
        x4^2 * x5 +
        x3 * x5^2 +
        x4 * x5^2 +
        x5^3 +
        x3^2 * x6 +
        x3 * x4 * x6 +
        x4^2 * x6 +
        x3 * x5 * x6 +
        x4 * x5 * x6 +
        x5^2 * x6 +
        x3 * x6^2 +
        x4 * x6^2 +
        x5 * x6^2 +
        x6^3,
        x4^4 +
        x4^3 * x5 +
        x4^2 * x5^2 +
        x4 * x5^3 +
        x5^4 +
        x4^3 * x6 +
        x4^2 * x5 * x6 +
        x4 * x5^2 * x6 +
        x5^3 * x6 +
        x4^2 * x6^2 +
        x4 * x5 * x6^2 +
        x5^2 * x6^2 +
        x4 * x6^3 +
        x5 * x6^3 +
        x6^4,
        x5^5 + x5^4 * x6 + x5^3 * x6^2 + x5^2 * x6^3 + x5 * x6^4 + x6^5,
        x6^6 - 1
    ]

    ku = Groebner.change_ordering(Groebner.ku10(ground=QQ), :degrevlex)
    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = gens(parent(first(ku)))
    gb = Groebner.groebner(ku, certify=true)

    @test gb == [
        x9 + 230965 // 116706 * x10 - 697789 // 116706,
        x8 + 92386 // 335511 * x10 + 578636 // 335511,
        x7 + 5173616 // 8563059 * x10 - 30862793 // 8563059,
        x6 + 87951472 // 192015891 * x10 + 488096201 // 192015891,
        x5 + 131980 // 18759 * x10 - 56944 // 18759,
        x4 + 32995 // 6084 * x10 - 63415 // 6084,
        x3 - 263960 // 229671 * x10 + 493631 // 229671,
        x2 - 65990 // 17043 * x10 + 151205 // 17043,
        x1 - 26396 // 3273 * x10 + 19850 // 3273,
        x10^2 - 71551 // 26396 * x10 + 45155 // 26396
    ]
end
