
using AbstractAlgebra

@testset "f4 adaptive coefficients" begin
    @test Groebner.determinechartype(2) == UInt8
    @test Groebner.determinechartype(100) == UInt16
    @test Groebner.determinechartype(2^15) == UInt32
    @test Groebner.determinechartype(2^31 - 1) == UInt64
    @test Groebner.determinechartype(2^60) == UInt128

    R, (x, y) = PolynomialRing(GF(2), ["x", "y"])
    fs = [
        x^2 + y
        x^3 + x * y
    ]
    @test Groebner.groebner(fs) == [x^2 + y]

    R, (x, y) = PolynomialRing(GF(37), ["x", "y"])
    fs = [
        x^2 - 4y
        x^3 - 3x * y
    ]
    @test Groebner.groebner(fs) == [y^2, x * y, x^2 + 33y]

    R, (x, y) = PolynomialRing(GF(1031), ["x", "y"])
    fs = [
        x^2 - 4y
        x^3 - 3x * y
    ]
    @test Groebner.groebner(fs) == [y^2, x * y, x^2 + 1027y]

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"])
    fs = [
        x^2 - 4y
        x^3 - 3x * y
    ]
    @test Groebner.groebner(fs) == [y^2, x * y, x^2 + 2147483643y]

    R, (x, y) = PolynomialRing(GF(2^40 + 15), ["x", "y"])
    fs = [
        x^2 - 4y
        x^3 - 3x * y
    ]
    @test Groebner.groebner(fs) == [y^2, x * y, x^2 + 1099511627787y]
end
