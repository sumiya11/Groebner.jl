
@testset "Normal form f4 finite fields" begin

    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"])

    # Regression test, Structural Identifiability
    @test Groebner.normalform([x], R(0)) == R(0)
    @test Groebner.normalform([x], R(1)) == R(1)

    Gs = [
        [x], [x, y], [x, y, z]
    ]

    @test Groebner.normalform(Gs[1], x) == 0
    @test Groebner.normalform(Gs[1], y) == y
    @test Groebner.normalform(Gs[3], x + y^2) == 0
    @test Groebner.normalform(Gs[2], 5x^2 - y^3 - 1) == -1
    @test Groebner.normalform(Gs[2], x + y + 3z) == 3z

    G = [
        x^2 + y*x + 1,
        y^2 - y
    ]

    @test Groebner.normalform(G, x + y + z) == x + y + z
    @test Groebner.normalform(G, x^2 + y^2) == y - y*x - 1
    @test Groebner.normalform(G, y^3) == Groebner.normalform(G, y^4) == y


    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x","y","z"], ordering=:degrevlex)
    G = [
        x^2 + y,
        y^2 + x
    ]
    @test Groebner.normalform(G, x^2 + y^2) == -x - y

    @test_nowarn Groebner.normalform([x, x*y], x)
end

@testset "Normal form -- a lot of variables" begin
    R, xs = PolynomialRing(QQ, ["x$i" for i in 1:63])

    GB = Groebner.groebner(xs)
    @test GB == reverse(xs)
    @test Groebner.normalform(GB, xs[1]) == R(0)


    R, xs = PolynomialRing(QQ, ["x$i" for i in 1:220])

    GB = Groebner.groebner(xs)
    @test GB == reverse(xs)
    @test Groebner.normalform(GB, xs[1]) == R(0)

    #=
    R, xs = PolynomialRing(QQ, ["x$i" for i in 1:1049])

    GB = Groebner.groebner(xs)
    @test GB == reverse(xs)
    @test Groebner.normalform(GB, xs[1]) == R(0)
    =#

end
