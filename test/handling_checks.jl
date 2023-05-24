
@testset "normalform checks" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"], ordering=:lex)

    @test Groebner.normalform([x, y], y, check=true) == R(0)
    @test_throws DomainError Groebner.normalform([x, x + 1], y, check=true)
    @test_throws DomainError Groebner.normalform([x, x + 1], [y], check=true)
end

@testset "fglm checks" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"], ordering=:lex)

    b = Groebner.fglm([x, y, z], check=true)
    @test b == reverse(collect(gens(parent(first(b)))))
    @test_throws DomainError Groebner.fglm([x, x + 1], check=true)
    @test_throws DomainError Groebner.fglm([x, x + 1], check=true)
end

@testset "kbase checks" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"], ordering=:lex)

    b = Groebner.kbase([x, y, z], check=true)
    @test b == [R(1)]
    @test_throws DomainError Groebner.kbase([x, x + 1], check=true)
    @test_throws DomainError Groebner.kbase([x, x + 1], check=true)
end
