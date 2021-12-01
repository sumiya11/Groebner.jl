

using .GroebnerBases: f4, rootn, katsura6

@testset "F4 over Finite Fields" begin

    ground = GF(2^31-1)
    R, (x, y, z, w) = PolynomialRing(ground, ["x", "y", "z", "w"])

    #########
    fs = [R(0)]
    Gzero = f4(fs)
    Gempty = f4(elem_type(R)[])

    # Root n
    ground = GF(2^31-1)
    fs = rootn(3, ground=ground)
    x1, x2, x3 = gens(parent(first(fs)))
    G = f4(fs, reduced=true)

    @test G == [x3^3 - 1, x2^2 + x2*x3 + x3^2, x1 + x2 + x3]

    #########
    fs = [
        x + y,
        x^2 + z^2,
        y^2 + w^2,
        x + 2*y + 3*z + 4*w
    ]
    G = f4(fs, reduced=true)

    @test G == [w^3, z*w + 894784854*w^2, z^2 - w^2,
                y + 3*z + 4*w, x - 3*z - 4*w]

    #########
    fs = [
        x*y^4 + y^3 + z*w^2,
        x^2 + y^4 + w^6,
        x*y*z*w
    ]
    G = f4(fs, reduced=true)
    @test G == [z^6*w^5 + z^2*w^15, y*z*w^11 + z^4*w^5,
                y*z^3*w^5 - z^2*w^9, y^2*z*w^7 + z^3*w^5,
                y^2*z^2*w^3 - y*z*w^7, y^3*z*w + z^2*w^3,
                y^9 + y^5*w^6 + y^3 + z*w^2, x*z^2*w^3,
                x*y*z*w, x*y^3 + x*z*w^2 - y^8 - y^4*w^6,
                x^2 + y^4 + w^6]

    

end
