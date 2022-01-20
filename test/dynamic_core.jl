
using DynamicPolynomials
@testset "DynamicPolynomials basic f4" begin
    @polyvar x y

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
