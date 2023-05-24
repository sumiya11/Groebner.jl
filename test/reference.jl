
using .Groebner:
    _normal_form_reference,
    _reducegb_reference,
    _reducegb_reference!,
    change_ordering,
    rootn

using AbstractAlgebra: ordering

@testset "reference normal form" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

    Gs = [[x], [x, y], [x, y, z]]

    @test _normal_form_reference(x, Gs[1]) == 0
    @test _normal_form_reference(y, Gs[1]) == y
    @test _normal_form_reference(x + y^2, Gs[3]) == 0
    @test _normal_form_reference(5x^2 - y^3 - 1, Gs[2]) == -1
    @test _normal_form_reference(x + y + 3z, Gs[2]) == 3z

    G = [x^2 + y * x + 1, y^2 - y]

    @test _normal_form_reference(x + y + z, G) == x + y + z
    @test _normal_form_reference(x^2 + y^2, G) == y - y * x - 1
    @test _normal_form_reference(y^3, G) == _normal_form_reference(y^4, G) == y

    G = [x * y + y, z + 5]
    @test _normal_form_reference(x * y, G) == -y
    @test _normal_form_reference(3z - 1, G) == -16

    @test _normal_form_reference(8x * y + 88z, [R(1)]) == 0

    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"], ordering=:degrevlex)
    G = [x^2 + y, y^2 + x]
    @test _normal_form_reference(x^2 + y^2, G) == -x - y
end

@testset "reference change order" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

    fs = [1x + 2y * z + 3z^3, y * z, 5x]
    fs_deg = change_ordering(fs, :degrevlex)
    R_deg = parent(first(fs_deg))

    @test ordering(R_deg) == :degrevlex
    @test collect(exponent_vectors(fs_deg[1])) == [[0, 0, 3], [0, 1, 1], [1, 0, 0]]
    @test collect(exponent_vectors(fs_deg[2])) == [[0, 1, 1]]
    @test collect(exponent_vectors(fs_deg[3])) == [[1, 0, 0]]

    fs_vv = change_ordering(fs, :lex)
    @test fs_vv == fs
end
