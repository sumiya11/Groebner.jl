
using DynamicPolynomials

@testset "Input-output DynamicPolynomials" begin

    @polyvar x y
    fs = [x, y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == [Polynomial(x), Polynomial(y)]

    fs = [x^2, 2y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == [Polynomial(x^2), Polynomial(2y)]

    fs = [(3//4)x, (2//9)y, -y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == [Polynomial((3//4)x), Polynomial((2//9)y), Polynomial(-y)]

    fs = [x^2*y + 3//4, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    fs = [x^2 + 1, (2^31 - 5)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs

    #= DynamicPolynomials does not work with AbstractAlgebra coefficients
    fs = [x^2*y + ff(1), ff(2)*x - (2^31 - 4)*y]
    ring, exps, cfs = Groebner.convert_to_internal(fs, :input)
    fsfs = Groebner.convert_to_output(ring, fs, exps, cfs)
    @test fsfs == fs
    =#

end

@testset "DynamicPolynomials io consistency" begin

end
