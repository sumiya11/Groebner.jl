
using .GroebnerBases: groebner


R, (x, y, z, w) = PolynomialRing(QQ, ["x", "y", "z", "w"])
Rz, (xz, yz, zz, wz) = PolynomialRing(ZZ, ["x", "y", "z", "w"])


@testset "Groebner with reconstruction" begin

    fs = [(12345678//12347)x,
          (222222221111123//2131232232097)y + z]
    G = GroebnerBases.groebner(fs)
    @test GroebnerBases.is_groebner(G, initial_gens=fs)


    fs = [
        (12345//12345678)x + y,
        x^2 + z^2,
        y^2 + w^2,
        x + 2*y + 3*z + 4*w
    ]
    G = GroebnerBases.groebner(fs)

    @test G == [w^3,
                z*w + 1692834523553//4063974*w^2,
                z^2 - 16935085031076//16933225*w^2,
                y - 12345//4106996*z - 4115//1026749*w,
                x + 6172839//2053498*z + 4115226//1026749*w]
    @test GroebnerBases.is_groebner(G, initial_gens=fs)

    # what if we are unlucky to start from initial prime in coefficients
    G = [
        (2^30 + 3)x,
        (1//(2^30 + 3))y
    ]
    G = GroebnerBases.groebner(G)
    @test G == [y, x]

end
