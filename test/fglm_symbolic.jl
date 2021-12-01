

using .GroebnerBases: compose_poly, decompose_poly, insert_nexts!

R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

@testset "fglm insert nexts" begin

    @test insert_nexts!(elem_type(R)[], x) == [z*x, y*x, x^2]
    @test insert_nexts!([x^3], x^2) == [x^2*z, x^2*y, x^3]
    @test insert_nexts!([x, y, z], R(1)) == [z, y, x]

end


@testset "fglm polys to vectors and vice versa" begin

    @test decompose_poly(-2x, [x]) == [QQ(-2)]

    f = R(1)
    monoms = [R(1), x, y, z]
    @test decompose_poly(f, monoms) == QQ.([1, 0, 0, 0])

    f = 2 + x + 3x^2 + -4z^2
    monoms = [R(1), x, x^2, y, y^2, z, z^2]
    @test decompose_poly(f, monoms) == QQ.([2, 1, 3, 0, 0, 0, -4])

    f = (2x + 3y + 4z)^3
    monoms = collect(monomials([f, (x^2 - 5y^3)^2, (1 - z)^2]))
    @test compose_poly( decompose_poly(f, monoms), monoms ) == f


end
