import AbstractAlgebra

@testset "output type inferred" begin
    R, (x, y) = polynomial_ring(AbstractAlgebra.QQ, ["x", "y"]; internal_ordering=:lex)
    @test @inferred Groebner.groebner([x, y]) == [y, x]
    @test @inferred Groebner.normalform([x, y], x + 1) == R(1)
    @test @inferred Groebner.normalform([x, y], [x + 1, y + 1, R(0)]) == [R(1), R(1), R(0)]
    @test @inferred Groebner.isgroebner([x, y]) == true

    R, (x, y) =
        polynomial_ring(AbstractAlgebra.GF(2^31 - 1), ["x", "y"]; internal_ordering=:degrevlex)
    @test @inferred Groebner.groebner([x, y]) == [y, x]
    @test @inferred Groebner.normalform([x, y], x + 1) == R(1)
    @test @inferred Groebner.normalform([x, y], [x + 1, y + 1, R(0)]) == [R(1), R(1), R(0)]
    @test @inferred Groebner.isgroebner([x, y]) == true

    context, gb = Groebner.groebner_learn([x, y])
    @test @inferred Groebner.groebner_learn([x, y])[2] == gb
    @test @inferred Groebner.groebner_apply!(context, [x, y]) == (true, gb)
end
