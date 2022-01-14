

@testset "Input-output AbstractAlgebra" begin
    using AbstractAlgebra

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs)
    fsfs = FastGroebner.export_basis(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs)
    fsfs = FastGroebner.export_basis(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs)
    fsfs = FastGroebner.export_basis(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs)
    fsfs = FastGroebner.export_basis(ring, fs, exps, cfs)
    @test fsfs == fs
end

@testset "Input-output generic :hasparent" begin
    using AbstractAlgebra

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs)
    ring.origring = :hasparent
    fsfs = FastGroebner.export_basis(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs)
    ring.origring = :hasparent
    fsfs = FastGroebner.export_basis(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs)
    ring.origring = :hasparent
    fsfs = FastGroebner.export_basis(ring, fs, exps, cfs)
    @test fsfs == fs

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:lex)
    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = FastGroebner.convert_to_internal(fs)
    ring.origring = :hasparent
    fsfs = FastGroebner.export_basis(ring, fs, exps, cfs)
    @test fsfs == fs
end

@testset "Input-output DynamicPolynomials" begin

end
