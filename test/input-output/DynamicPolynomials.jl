using DynamicPolynomials

@testset "io consistency dynamicpolynomials" begin
    @polyvar x y

    fs = [x, y]
    @test Groebner.groebner(fs) isa Vector{Polynomial{true, Int64}}

    fs = [UInt16(2)x + UInt16(3), UInt16(2)y]
    @test Groebner.groebner(fs) isa Vector{Polynomial{true, UInt16}}

    fs = [BigInt(1)x + BigInt(20), y]
    @test Groebner.groebner(fs) isa Vector{Polynomial{true, BigInt}}

    fs = [x + BigInt(20) // BigInt(3), y]
    @test Groebner.groebner(fs) isa Vector{Polynomial{true, Rational{BigInt}}}

    fs = [34343343433x * y^2 + 3431234567833, 3434343434x * y - 342343242342]
    @test_throws DomainError Groebner.groebner(fs)

    fs = [BigInt(34343343433)x * y^2 + 3431234567833, 3434343434x * y - 342343242342]
    @test Groebner.groebner(fs) isa Vector{Polynomial{true, BigInt}}
end
