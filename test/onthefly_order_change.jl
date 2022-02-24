
@testset "check parameter order change" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x","y","z"], ordering=:lex)

    @test Groebner.groebner([y^2, x], ordering=:lex) == [y^2, x]
    @test Groebner.groebner([y^2, x], ordering=:degrevlex) == [x, y^2]
    @test Groebner.groebner([y^2, x], ordering=:deglex) == [x, y^2]

    @test Groebner.groebner([y^2, x], forsolve=true) == [y^2, x]

end
