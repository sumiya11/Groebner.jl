
using .GroebnerBases: fglm, groebner, f4

@testset "FGLM over Finite Fields" begin

    ground = GF(2^31-1)
    R, (x, y, z) = PolynomialRing(ground, ["x", "y", "z"], ordering=:degrevlex)

    fs_deg = [
        x^2 + 2y^2 - y - 2z,
        x^2 - 8y^2 + 10z - 1,
        x^2 - 7y*z
    ]

    gb_deg = f4(fs_deg)
    println(gb_deg)

    gb_lex = fglm(gb_deg)
    println(gb_lex)
end
