

using AbstractAlgebra

@testset "Input-output AbstractAlgebra" begin

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:deglex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:deglex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:deglex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:deglex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs
end

@testset "Input-output generic :hasparent" begin

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    ring.origring = :hasparent
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    ring.origring = :hasparent
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    ring.origring = :hasparent
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    ring.origring = :hasparent
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs
end

# using DynamicPolynomials

@testset "Input-output DynamicPolynomials" begin

    #=
    ff = GF(2^31 - 1)
    @polyvar x y

    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs
    =#
    
    #=
    fs = [x^2 + 1, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs
    =#

    #= DynamicPolynomials does not work with AbstractAlgebra
    fs = [x^2*y + ff(1), ff(2)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs, :input)
    fsfs = FastGroebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs
    =#

end
