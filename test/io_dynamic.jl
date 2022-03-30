
using DynamicPolynomials

@testset "Input-output DynamicPolynomials" begin

    @polyvar x y
    fs = [x, y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    meta = Groebner.set_metaparameters(ring, :input, false, false, :exact)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
    @test fsfs == [Polynomial(x), Polynomial(y)]

    fs = [x^2, 2y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    meta = Groebner.set_metaparameters(ring, :input, false, false, :exact)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
    @test fsfs == [Polynomial(x^2), Polynomial(2y)]

    fs = [(3//4)x, (2//9)y, -y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    meta = Groebner.set_metaparameters(ring, :input, false, false, :exact)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
    @test fsfs == [Polynomial((3//4)x), Polynomial((2//9)y), Polynomial(-y)]

    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    meta = Groebner.set_metaparameters(ring, :input, false, false, :exact)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
    @test fsfs == fs

    fs = [x^2 + 1, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    meta = Groebner.set_metaparameters(ring, :input, false, false, :exact)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
    @test fsfs == fs

    #= DynamicPolynomials does not work with AbstractAlgebra coefficients
    fs = [x^2*y + ff(1), ff(2)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs
    =#

end

@testset "DynamicPolynomials io consistency" begin
    @polyvar x y

    fs = [x, y]
    @test Groebner.groebner(fs) isa Vector{Polynomial{true, Int64}}

    fs = [UInt16(2)x + UInt16(3), UInt16(2)y]
    @test Groebner.groebner(fs) isa Vector{Polynomial{true, UInt16}}

    fs = [BigInt(1)x + BigInt(20), y]
    @test Groebner.groebner(fs) isa Vector{Polynomial{true, BigInt}}

    fs = [x + BigInt(20)//BigInt(3), y]
    @test Groebner.groebner(fs) isa Vector{Polynomial{true, Rational{BigInt}}}

    fs = [34343343433x*y^2 + 3431234567833, 3434343434x*y - 342343242342]
    @test_throws DomainError Groebner.groebner(fs)

    fs = [BigInt(34343343433)x*y^2 + 3431234567833, 3434343434x*y - 342343242342]
    @test Groebner.groebner(fs) isa Vector{Polynomial{true, BigInt}}
end
