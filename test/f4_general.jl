

using .GroebnerBases: f4, rootn, katsura6, is_groebner


function f4_tests(; f4_config...)
    ground = GF(2^31 - 1)
    R, (x, y, z, w) = PolynomialRing(ground, ["x", "y", "z", "w"])

    #########
    # I am not entirely sure what do we expect here
    Gzero = f4(zeros(R, 1); f4_config...)
    Gempty = f4(zeros(R, 0); f4_config...)

    # Root n
    ground = GF(2^31-1)
    fs = rootn(3, ground=ground)
    x1, x2, x3 = gens(parent(first(fs)))
    G = f4(fs; f4_config...)

    @test G == [x3^3 - 1, x2^2 + x2*x3 + x3^2, x1 + x2 + x3]

    @test is_groebner(G, initial_gens=fs)
    @test !is_groebner( append!(G, [x3]) , initial_gens=fs)

    G = [x3^3 - 1, x2^2 + x2*x3 + x3^2, x1 + x2 + x3]
    # we want to omit duplicates for now, broken test !
    # upd: the test is okay
    @test is_groebner( append!(G, [x1 + x2 + x3]) , initial_gens=fs)

    #########
    fs = [
        x + y,
        x^2 + z^2,
        y^2 + w^2,
        x + 2*y + 3*z + 4*w
    ]
    G = f4(fs; f4_config...)

    @test G == [w^3, z*w + 894784854*w^2, z^2 - w^2,
                y + 3*z + 4*w, x - 3*z - 4*w]
    @test is_groebner(G, initial_gens=fs)

    #########
    fs = [
        x*y^4 + y^3 + z*w^2,
        x^2 + y^4 + w^6,
        x*y*z*w
    ]
    G = f4(fs; f4_config...)
    @test G == [z^6*w^5 + z^2*w^15, y*z*w^11 + z^4*w^5,
                y*z^3*w^5 - z^2*w^9, y^2*z*w^7 + z^3*w^5,
                y^2*z^2*w^3 - y*z*w^7, y^3*z*w + z^2*w^3,
                y^9 + y^5*w^6 + y^3 + z*w^2, x*z^2*w^3,
                x*y*z*w, x*y^3 + x*z*w^2 - y^8 - y^4*w^6,
                x^2 + y^4 + w^6]
    @test is_groebner(G)
end


@testset ":dense F4 over FF in :lex" begin

    f4_tests(; linalg=:dense)
end

@testset ":sparse F4 over FF in :lex" begin

    f4_tests(; linalg=:sparse)
end
