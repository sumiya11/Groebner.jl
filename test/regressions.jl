
@testset "monomial overflow" begin
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"], ordering=:degrevlex)
    for monoms in [:default, :packed, :sparse]
        gb_1 = [x * y^100 + y, x^100 * y + y^100, y^199 + 2147483646 * x^99 * y]
        gb_2 = [x * y^200 + y, x^200 * y + y^200, y^399 + 2147483646 * x^199 * y]
        gb_3 = [x * y^1000 + y, x^1000 * y + y^1000, y^1999 + 2147483646 * x^999 * y]
        @test Groebner.groebner([x^100 * y + y^100, x * y^100 + y], monoms=monoms) == gb_1
        @test Groebner.groebner([x^200 * y + y^200, x * y^200 + y], monoms=monoms) == gb_2
        @test Groebner.groebner([x^1000 * y + y^1000, x * y^1000 + y], monoms=monoms) ==
              gb_3

        @test Groebner.isgroebner(gb_1)
        @test Groebner.isgroebner(gb_2)
        @test Groebner.isgroebner(gb_3)

        @test Groebner.normalform(gb_1, [x, y, R(1), R(0), x^1000]) ==
              [x, y, R(1), R(0), x^1000]
        @test Groebner.normalform(gb_2, [x, y, R(1), R(0), x^1000]) ==
              [x, y, R(1), R(0), x^1000]
        @test Groebner.normalform(gb_3, [x, y, R(1), R(0), x^10]) ==
              [x, y, R(1), R(0), x^10]
        @test Groebner.normalform(gb_3, [x, y, R(1), R(0), x^1000]) ==
              [x, y, R(1), R(0), x^1000]
    end
end

@testset "SI.jl normalform bug" begin
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"])

    @test Groebner.normalform([x], R(0)) == R(0)
    @test Groebner.normalform([x], R(1)) == R(1)
    @test Groebner.normalform([x], [R(0)]) == [R(0)]
    @test Groebner.normalform([x], [R(1)]) == [R(1)]
    @test Groebner.normalform([x], [R(0), R(1), R(0)]) == [R(0), R(1), R(0)]
end

@testset "SI.jl cmp bug" begin
    # this may crash if the comparator is invalid
    function bug(gens::Vector{Int}, exps::Vector{Vector{T}}) where {T}
        inds = collect(1:length(gens))
        cmp  = (x, y) -> Groebner.monom_isless(exps[gens[x]], exps[gens[y]], Groebner.DegRevLex())
        sort!(inds, lt=cmp)
        @test true
    end

    gens = [
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        3,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        3,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        3,
        2,
        2,
        2,
        2,
        2
    ]

    exps = Groebner.ExponentVector{UInt64}[
        [0x0000, 0x0000],
        [0x0002, 0x0002],
        [0x0002, 0x0002]
    ]

    bug(gens, exps)
end
