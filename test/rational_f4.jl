
using .FastGroebner: groebner


@testset "Groebner modular" begin
    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)

    fs = [x, y]
    G = FastGroebner.groebner(fs)
    @test G ≂ [y, x]

    fs = [5y, 3//4*x]
    G = FastGroebner.groebner(fs)
    @test G ≂ [y, x]

    fs = [x^2 - (1//6)*y, x*y]
    G = FastGroebner.groebner(fs)
    @test G ≂ [y^2, x*y, x^2 - (1//6)y]

    fs = [QQ(11, 3)*x^2*y - QQ(2, 4)*x*y - QQ(1, 7)*y,
            QQ(1, 1)*x*y + QQ(7, 13)*x]
    G = FastGroebner.groebner(fs)
    @test G ≂ [y^2 + 7//13*y,
                x*y + 7//13*x,
                x^2 - 3//22*x + 39//539*y]

end

@testset "Groebner modular corner cases" begin
    R, (x, y, z, w) = PolynomialRing(QQ, ["x", "y", "z", "w"], ordering=:degrevlex)

    fs = [(12345678//12347)x,
          (222222221111123//2131232232097)y + z]
    G = FastGroebner.groebner(fs)
    @test G ≂ [y + 2131232232097//222222221111123*z, x]

    fs = [
        (12345//12345678)x + y,
        x^2 + z^2,
        y^2 + w^2,
        x + 2*y + 3*z + 4*w
    ]
    G = FastGroebner.groebner(fs)

    @test G ≂ [w^3,
                z*w + 1692834523553//4063974*w^2,
                z^2 - 16935085031076//16933225*w^2,
                y - 12345//4106996*z - 4115//1026749*w,
                x + 6172839//2053498*z + 4115226//1026749*w]
    @test FastGroebner.isgroebner(G, initial_gens=fs)

    # what if we are unlucky to start from initial prime in coefficients
    G = [
        (2^30 + 3)x,
        (1//(2^30 + 3))y
    ]
    G = FastGroebner.groebner(G)
    @test G ≂ [y, x]

    G = [
        (2^31 - 1)*x + y,
        (1//(2^31 - 1))*y + 1
    ]
    G = FastGroebner.groebner(G)
    @test G ≂ [y - 12, x - 1]

end
