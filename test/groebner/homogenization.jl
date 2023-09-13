using AbstractAlgebra
using Combinatorics
import Random

@testset "homogenization simple" begin
    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)

    @test Groebner.groebner([R(3)], homogenize=:yes) == [R(1)]
    @test Groebner.groebner([R(0)], homogenize=:yes) == [R(0)]
    @test Groebner.groebner([x], homogenize=:yes) == [x]
    @test Groebner.groebner([x + 11y], homogenize=:yes) == [x + 11y]
    @test Groebner.groebner([x, y], homogenize=:yes) == [y, x]
    @test Groebner.groebner([x + 1, y + 2], homogenize=:yes) == [y + 2, x + 1]

    # some tests with homogenize=:yes


    # compare all three values for homogenize
    for case in []
        
    end
end