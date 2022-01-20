
using AbstractAlgebra

@testset "Input-output AbstractAlgebra" begin

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:deglex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:deglex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:deglex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:deglex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y, x*y - y^2, x*y - x^2]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs
end

@testset "Input-output generic :hasparent" begin

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    ring.origring = :hasparent
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    ring.origring = :hasparent
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    ring.origring = :hasparent
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    ring.origring = :hasparent
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs
end
