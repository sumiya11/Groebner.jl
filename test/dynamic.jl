
using DynamicPolynomials
@testset "DynamicPolynomials basic groebner" begin
    @polyvar x y

    @test Groebner.groebner([x, y]) == Groebner.groebner([y, x]) == [y, x]
    @test Groebner.groebner([5x, y]) == Groebner.groebner([y, 5x]) == [y, x]
    @test Groebner.groebner([x, y, y]) == Groebner.groebner([y, y, x]) == [y, x]
    @test Groebner.groebner([x + 3//2, y]) == [y, x + 3//2]
    @test Groebner.groebner([x + 3, 2y - 6]) == [y - 3, x + 3]

    fs = [x^2 - 5y]
    G = Groebner.groebner(fs)
    @test G == [x^2 - 5y]

    fs = [
        x + y^2,
        x*y - y^2
    ]
    G = Groebner.groebner(fs)
    @test Groebner.isgroebner(G)

end
