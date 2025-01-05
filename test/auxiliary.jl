using AbstractAlgebra
using Test

@testset "lead" begin
    R, (x, y) = QQ["x", "y"]

    @test Groebner.lead(2x - 3y) == 2x
    @test Groebner.lead(2x - 3y, ordering=Lex(x, y)) == 2x
    @test Groebner.lead(2x - 3y, ordering=Lex(y, x)) == -3y

    R, (X, Y) = fraction_field(R)["X", "Y"]
    @test Groebner.lead(2 * x * Y - 3 * y * X, ordering=Lex(X, Y)) == -3 * y * X
    @test Groebner.lead(2 * x * Y - 3 * y * X, ordering=Lex(Y, X)) == 2 * x * Y
    @test Groebner.lead(zero(R)) == zero(R)
    @test Groebner.lead(-one(R)) == -one(R)
end
