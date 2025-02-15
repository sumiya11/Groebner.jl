using AbstractAlgebra
using Test

@testset "leading term" begin
    R, (x, y) = QQ["x", "y"]

    @test Groebner.leading_term(2x - 3y) == 2x
    @test Groebner.leading_term(2x - 3y, ordering=Lex(x, y)) == 2x
    @test Groebner.leading_term(2x - 3y, ordering=Lex(y, x)) == -3y

    @test_throws DomainError Groebner.leading_term([1, 2, 3])

    R, (X, Y) = fraction_field(R)["X", "Y"]
    @test Groebner.leading_term(2 * x * Y - 3 * y * X, ordering=Lex(X, Y)) == -3 * y * X
    @test Groebner.leading_term(2 * x * Y - 3 * y * X, ordering=Lex(Y, X)) == 2 * x * Y
    @test Groebner.leading_term(zero(R)) == zero(R)
    @test Groebner.leading_term(-one(R)) == -one(R)
end

@testset "leading ideal" begin
    R, (x, y) = QQ["x", "y"]
    @test Groebner.leading_ideal([R(5)]) == [R(1)]
    @test Groebner.leading_ideal([x - 5, y - 1]) == [y, x]
    @test Groebner.leading_ideal([x * y, y, x, R(0)]) == [y, x]
    @test Groebner.leading_ideal([x^1000]) == [x^1000]
    @test Groebner.leading_ideal([R(0), x, R(0)]) == [x]

    @test_throws DomainError Groebner.leading_ideal([1, 2, 3])

    x0, x1, x2, x3 = gens(parent(Groebner.Examples.katsuran(3)[1]))
    @test Groebner.leading_ideal(Groebner.Examples.katsuran(3), ordering=Lex()) ==
          [x3^8, x2, x1, x0]

    x0, x1, x2, x3 = gens(parent(Groebner.Examples.katsuran(3)[1]))
    @test Groebner.leading_ideal(Groebner.Examples.katsuran(3), ordering=DegRevLex()) ==
          [x0, x2^2, x1 * x2, x1^2, x2 * x3^2, x1 * x3^2, x3^4]
end

@testset "quotient basis" begin
    R, (x, y) = QQ["x", "y"]
    @test Groebner.quotient_basis([R(5)]) == Vector{elem_type(R)}()
    @test Groebner.quotient_basis([x - 5, y - 1]) == [R(1)]
    @test Groebner.quotient_basis([x * y, y, x, R(0)]) == [R(1)]
    @test length(Groebner.quotient_basis([x^100, y^2])) == 100 * 2
    @test Groebner.quotient_basis([R(0), x, R(0), y]) == [R(1)]

    @test_throws DomainError Groebner.quotient_basis([])
    @test_throws DomainError Groebner.quotient_basis([x])
    @test_throws DomainError Groebner.quotient_basis([R(0)])

    for k in [
        AbstractAlgebra.GF(2^30 + 3),
        fraction_field(AbstractAlgebra.QQ["t"][1]),
        AbstractAlgebra.QQ
    ]
        for (i, sys) in enumerate([
            Groebner.Examples.katsuran(1, k=k),
            Groebner.Examples.katsuran(2, k=k),
            Groebner.Examples.katsuran(3, k=k)
        ])
            @test length(Groebner.quotient_basis(sys)) == 2^i
        end
    end
end

@testset "dimension" begin
    R, (x, y) = QQ["x", "y"]
    @test Groebner.dimension([R(5)]) == -1
    @test Groebner.dimension([x - 5, y - 1]) == 0
    @test Groebner.dimension([x - 5]) == 1
    @test Groebner.dimension([R(0)]) == -1
    @test Groebner.dimension([x * y, y, x, R(0)]) == 0

    @test_throws DomainError Groebner.dimension([])

    n = 100
    R, x = polynomial_ring(AbstractAlgebra.GF(2^30 + 3), ["x$i" for i in 1:n])
    @test Groebner.dimension([sum(x)]) == n - 1
    @test Groebner.dimension([sum(x), prod(x), sum([i for i in 1:n] .* x)]) == n - 3
end
