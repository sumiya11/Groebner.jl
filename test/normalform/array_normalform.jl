
using AbstractAlgebra

@testset "normalform of an array" begin
    for field in [GF(17), GF(2^31 - 1), ZZ, QQ]
        R, (x, y) = PolynomialRing(field, ["x", "y"])
        gb = [x, y]
        @test Groebner.normalform(gb, [x, y + 1]) == [R(0), R(1)]
        @test Groebner.normalform(gb, [y + 1, x]) == [R(1), R(0)]
        @test Groebner.normalform(
            gb,
            [x, 3y, 4x + 5y, 6x * y, 7x + 1, y + 2, x^2 + y^2 + 3]
        ) == [R(0), R(0), R(0), R(0), R(1), R(2), R(3)]

        gb = [x]
        @test Groebner.normalform(
            gb,
            [x, 3y, 4x + 5y, 6x * y, 7x + 1, y + 2, x^2 + y^2 + 3]
        ) == [R(0), 3y, 5y, R(0), R(1), y + 2, y^2 + 3]

        gb = [x^2 + y^2, y^3 - y]
        @test Groebner.normalform(
            gb,
            [R(1), R(5), R(16), 3x, y, 4x * y + 2, y^2, x * y^2 - 8x + y]
        ) == [R(1), R(5), R(16), 3x, y, 4x * y + 2, y^2, x * y^2 - 8x + y]
    end
end
