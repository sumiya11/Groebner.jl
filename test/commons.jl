
using .Groebner: normal_form, reducegb, reducegb!,
                    change_ordering, rootn

@testset "Normal Form" begin

    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

    Gs = [
        [x], [x, y], [x, y, z]
    ]

    @test normal_form(x, Gs[1]) == 0
    @test normal_form(y, Gs[1]) == y
    @test normal_form(x + y^2, Gs[3]) == 0
    @test normal_form(5x^2 - y^3 - 1, Gs[2]) == -1
    @test normal_form(x + y + 3z, Gs[2]) == 3z

    G = [
        x^2 + y*x + 1,
        y^2 - y
    ]

    @test normal_form(x + y + z, G) == x + y + z
    @test normal_form(x^2 + y^2, G) == y - y*x - 1
    @test normal_form(y^3, G) == normal_form(y^4, G) == y

    G = [
        3x*y + y,
        4z + 5,
    ]
    @test normal_form(x*y, G) == -y//3
    @test normal_form(3z - 1, G) == -19//4

    @test normal_form(8x*y + 88z, [R(5)]) == 0

    R, (x, y, z) = PolynomialRing(QQ, ["x","y","z"], ordering=:degrevlex)
    G = [
        x^2 + y,
        y^2 + x
    ]
    @test normal_form(x^2 + y^2, G) == -x - y

end

@testset "Change ordering" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

    fs = [1x + 2y*z + 3z^3, y*z, 5x]
    fs_deg = change_ordering(fs, :degrevlex)
    R_deg = parent(first(fs_deg))

    @test ordering(R_deg) == :degrevlex
    @test collect(exponent_vectors(fs_deg[1])) == [[0, 0, 3], [0, 1, 1], [1, 0, 0]]
    @test collect(exponent_vectors(fs_deg[2])) == [[0, 1, 1]]
    @test collect(exponent_vectors(fs_deg[3])) == [[1, 0, 0]]

    fs_vv = change_ordering(fs, :lex)
    @test fs_vv == fs
end
