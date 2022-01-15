
#=
    Here we check that the result of computation is
    of the desired type in respect to input types
=#

using AbstractAlgebra

@testset "AbstractAlgebra io consistency" begin
    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    gb = FastGroebner.groebner(fs)
    @test parent(gb[1]) == R

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3//2, (2^31 - 5)*x - (2^31 - 4)*y]
    gb = FastGroebner.groebner(fs)
    @test parent(gb[1]) == R

    R, (x, y) = PolynomialRing(ZZ, ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    gb = FastGroebner.groebner(fs)
    @test parent(gb[1]) == R

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    gb = FastGroebner.groebner(fs)
    @test parent(gb[1]) == R

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3//2, (2^31 - 5)*x - (2^31 - 4)*y]
    gb = FastGroebner.groebner(fs)
    @test parent(gb[1]) == R

    R, (x, y) = PolynomialRing(ZZ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    gb = FastGroebner.groebner(fs)
    @test parent(gb[1]) == R
end

@testset "MultivariatePolynomials io consistency" begin

end
