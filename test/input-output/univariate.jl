
using AbstractAlgebra

@testset "aa univariate io" begin
    R, x = PolynomialRing(GF(2^31 - 1), "x")
    @test Groebner.groebner([x^2 - 4, x + 2]) == [x + 2]

    R, x = PolynomialRing(QQ, "x")
    @test Groebner.groebner([x^2 - 4, x + 2]) == [x + 2]

    R, x = PolynomialRing(ZZ, "x")
    @test Groebner.groebner([x^2 - 4, x + 2]) == [x + 2]
end
