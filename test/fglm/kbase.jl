
@testset "kbase simple" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering=:degrevlex)

    @test Groebner.kbase([x, y]) == [R(1)]
    @test Groebner.kbase([x + y, y]) == [R(1)]
    @test Groebner.kbase([-2x + 8, 10y + 1]) == [R(1)]

    @test Groebner.kbase([x^2 + 1, y^2 - x - 1]) == [R(1), y, x, x * y]
    @test Groebner.kbase([8x^5 + 1, y - 1]) == [R(1), x, x^2, x^3, x^4]
    @test Groebner.kbase([x^2, y^2, x * y]) == [R(1), y, x]
    @test Groebner.kbase([x^3, y^3, x * y]) == [R(1), y, x, y^2, x^2]
    @test length(Groebner.kbase([8x^77 + 1, y^34 - 1])) == 77 * 34

    R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering=:lex)
    @test Groebner.kbase([x, y]) == [R(1)]
    @test Groebner.kbase([x + y, y]) == [R(1)]
    @test Groebner.kbase([-2x + 8, 10y + 1]) == [R(1)]

    @test Groebner.kbase([x^2 + 1, y^2 - 1]) == [R(1), y, x, x * y]
    @test Groebner.kbase([8x^5 + 1, y - 1]) == [R(1), x, x^2, x^3, x^4]
    @test Groebner.kbase([x^2, y^2, x * y]) == [R(1), y, x]
    @test Groebner.kbase([x^3, y^3, x * y]) == [R(1), y, y^2, x, x^2]
    @test length(Groebner.kbase([8x^77 + 1, y^34 - 1])) == 77 * 34

    R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering=:deglex)

    @test Groebner.kbase([x, y]) == [R(1)]
    @test Groebner.kbase([x + y, y]) == [R(1)]
    @test Groebner.kbase([-2x + 8, 10y + 1]) == [R(1)]

    @test Groebner.kbase([x^2 + 1, y^2 - x - 1]) == [R(1), y, x, x * y]
    @test Groebner.kbase([8x^5 + 1, y - 1]) == [R(1), x, x^2, x^3, x^4]
    @test Groebner.kbase([x^2, y^2, x * y]) == [R(1), y, x]
    @test Groebner.kbase([x^3, y^3, x * y]) == [R(1), y, x, y^2, x^2]
    @test length(Groebner.kbase([8x^77 + 1, y^34 - 1])) == 77 * 34
end

@testset "kbase zeros" begin
    for field in [GF(2), GF(2^31 - 1), QQ]
        R, (x, y, z) = polynomial_ring(field, ["x", "y", "z"], ordering=:lex)

        @test_throws DomainError Groebner.kbase([R(0)])
        @test_throws DomainError Groebner.kbase([R(0), R(0), R(0)])

        @test Groebner.kbase([x, y^2, R(0), z, R(0)]) == [R(1), y]
    end
end

@testset "kbase big" begin
    root4 = Groebner.rootn(4)
    gb = Groebner.groebner(root4)
    R = parent(first(root4))
    (x1, x2, x3, x4) = gens(R)
    @test length(Groebner.kbase(gb)) == 24

    root6 = Groebner.rootn(6)
    gb = Groebner.groebner(root6)
    @test length(Groebner.kbase(gb)) == 720
    root6 = Groebner.rootn(6)
    gb = Groebner.groebner(root6, ordering=Groebner.DegRevLex())
    @test length(Groebner.kbase(gb)) == 720

    noon2 = Groebner.noonn(2, ordering=:degrevlex)
    R = parent(first(noon2))
    (x1, x2) = gens(R)
    gb = Groebner.groebner(noon2, ordering=Groebner.DegRevLex())
    @test Groebner.kbase(gb) == [R(1), x2, x1, x2^2, x1 * x2]

    noon = Groebner.noonn(3)
    R = parent(first(noon))
    gb = Groebner.groebner(noon, ordering=Groebner.DegRevLex())
    @test length(Groebner.kbase(gb, ordering=Groebner.DegRevLex())) == 21

    sys = Groebner.eco10(k=GF(2^31 - 1), ordering=:degrevlex)
    R = parent(first(sys))
    gb = Groebner.groebner(sys, ordering=Groebner.DegRevLex())
    @test length(Groebner.kbase(gb, ordering=Groebner.DegRevLex())) == 256
end
