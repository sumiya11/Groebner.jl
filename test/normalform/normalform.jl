
@testset "normalform" begin
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"])

    # Regression test, Structural Identifiability
    @test Groebner.normalform([x], R(0)) == R(0)
    @test Groebner.normalform([x], R(1)) == R(1)

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

    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"], ordering=:degrevlex)
    G = [x^2 + y, y^2 + x]
    @test Groebner.normalform(G, x^2 + y^2) == -x - y

    @test_nowarn Groebner.normalform([x, x * y], x)
end

@testset "normalform many variables" begin
    R, xs = PolynomialRing(QQ, ["x$i" for i in 1:63])

    GB = Groebner.groebner(xs)
    @test GB == reverse(xs)
    @test Groebner.normalform(GB, xs[1]) == R(0)

    R, xs = PolynomialRing(QQ, ["x$i" for i in 1:220])

    GB = Groebner.groebner(xs)
    @test GB == reverse(xs)
    @test Groebner.normalform(GB, xs[1]) == R(0)
end

@testset "normalform of an array" begin
    for field in [GF(17), GF(2^31 - 1), ZZ, QQ]
        R, (x, y) = PolynomialRing(field, ["x", "y"])
        gb = [x, y]
        @test Groebner.normalform(gb, [x, y + 1]) == [R(0), R(1)]
        @test Groebner.normalform(gb, [y + 1, x]) == [R(1), R(0)]
        @test Groebner.normalform(
            gb,
            [x, 3y, 4x + 5y, 6x * y, 7x + 1, y + 2, x^2 + y^2 + 3]
        ) == [R(0), R(0), R(0), R(0), R(1), R(2), R(3)]

        gb = [x]
        @test Groebner.normalform(
            gb,
            [x, 3y, 4x + 5y, 6x * y, 7x + 1, y + 2, x^2 + y^2 + 3]
        ) == [R(0), 3y, 5y, R(0), R(1), y + 2, y^2 + 3]

        gb = [x^2 + y^2, y^3 - y]
        @test Groebner.normalform(
            gb,
            [R(1), R(5), R(16), 3x, y, 4x * y + 2, y^2, x * y^2 - 8x + y]
        ) == [R(1), R(5), R(16), 3x, y, 4x * y + 2, y^2, x * y^2 - 8x + y]
    end
end
