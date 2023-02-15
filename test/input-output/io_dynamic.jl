
using DynamicPolynomials
using Random

@testset "input-output dynamicpolynomials" begin
    representation = Groebner.default_safe_representation(Groebner.NotPacked{UInt64}())
    rng = Random.MersenneTwister(42)
    for representation in [
            Groebner.default_safe_representation(Groebner.NotPacked{UInt64}()),
            Groebner.Representation{Groebner.PackedPair2{UInt64, UInt16}}()
        ]

        @polyvar x y
        fs = [x, y]
        ring, exps, cfs = Groebner.convert_to_internal(representation, fs, Groebner.InputOrdering())
        meta = Groebner.set_metaparameters(ring, Groebner.InputOrdering(), false, false, :exact, rng)
        fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
        @test fsfs == [Polynomial(x), Polynomial(y)]

        fs = [x^2, 2y]
        ring, exps, cfs = Groebner.convert_to_internal(representation, fs, Groebner.InputOrdering())
        meta = Groebner.set_metaparameters(ring, Groebner.InputOrdering(), false, false, :exact, rng)
        fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
        @test fsfs == [Polynomial(x^2), Polynomial(2y)]

        fs = [(3//4)x, (2//9)y, -y]
        ring, exps, cfs = Groebner.convert_to_internal(representation, fs, Groebner.InputOrdering())
        meta = Groebner.set_metaparameters(ring, Groebner.InputOrdering(), false, false, :exact, rng)
        fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
        @test fsfs == [Polynomial((3//4)x), Polynomial((2//9)y), Polynomial(-y)]

        fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
        ring, exps, cfs = Groebner.convert_to_internal(representation, fs, Groebner.InputOrdering())
        meta = Groebner.set_metaparameters(ring, Groebner.InputOrdering(), false, false, :exact, rng)
        fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
        @test fsfs == fs

        fs = [x^2 + 1, (2^31 - 5)*x - (2^31 - 4)*y]
        ring, exps, cfs = Groebner.convert_to_internal(representation, fs, Groebner.InputOrdering())
        meta = Groebner.set_metaparameters(ring, Groebner.InputOrdering(), false, false, :exact, rng)
        fsfs = Groebner.convert_to_output(ring, fs, exps, cfs, meta)
        @test fsfs == fs

        #= DynamicPolynomials does not work with AbstractAlgebra coefficients
        fs = [x^2*y + ff(1), ff(2)*x - (2^31 - 4)*y]
        ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
        fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
        @test fsfs == fs
        =#
    end
end

@testset "io consistency dynamicpolynomials" begin
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
