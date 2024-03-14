
@testset "groebner, change matrix" begin
    R, (x, y, z) =
        polynomial_ring(GF(2^30 + 3), ["x", "y", "z"], internal_ordering=:degrevlex)
    f = [x * y * z - 1, x * y + x * z + y * z, x + y + z]
    g, m = Groebner.groebner_with_change_matrix(f)
    @test m * f == g
    @test size(m) == (length(g), length(f))
    @test Groebner.isgroebner(g)

    for ground in [GF(2^30 + 3), QQ]
        R, (x, y) = polynomial_ring(ground, ["x", "y"], internal_ordering=:degrevlex)

        cases = [
            [R(0), R(0), R(0)],
            [R(0), R(1), R(0)],
            [x],
            [y],
            [x, x, x, x, R(0), R(0), x, x, x, x],
            [x + y, R(1), x^5 - y^5],
            [x^200 + y^250, x^100 * y^200 + 1],
            Groebner.katsuran(4, k=ground),
            Groebner.noonn(4, k=ground),
            Groebner.cyclicn(4, k=ground),
            [x^i * y + x for i in 1:500]
        ]

        for f in cases
            g, m = Groebner.groebner_with_change_matrix(f)
            @test m * f == g
            @test size(m) == (length(g), length(f))
            @test Groebner.isgroebner(g)
        end
    end

    f = Groebner.katsuran(4, k=QQ)
    g, m = Groebner.groebner_with_change_matrix(f, ordering=Groebner.DegRevLex())
    @test m * f == g
    @test_throws DomainError Groebner.groebner_with_change_matrix(
        f,
        ordering=Groebner.DegLex()
    )
    @test_throws DomainError Groebner.groebner_with_change_matrix(f)
end
