
@testset "regression, SI.jl normalform" begin
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"])

    @test Groebner.normalform([x], R(0)) == R(0)
    @test Groebner.normalform([x], R(1)) == R(1)
    @test Groebner.normalform([x], [R(0)]) == [R(0)]
    @test Groebner.normalform([x], [R(1)]) == [R(1)]
    @test Groebner.normalform([x], [R(0), R(1), R(0)]) == [R(0), R(1), R(0)]
end

@testset "regression, ordering of empty" begin
    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)

    ord = Groebner.DegRevLex()
    gb1 = Groebner.groebner([x, y], ordering=ord)
    gb2 = Groebner.groebner([R(0)], ordering=ord)
    @test parent(first(gb1)) == R
    @test parent(first(gb2)) == R
    nf1 = Groebner.normalform([x, y], x, ordering=ord)
    nf2 = Groebner.normalform([x, y], R(0), ordering=ord)
    @test parent(nf1) == R
    @test parent(nf2) == R
end

@testset "regression, SI.jl cmp" begin
    # this may crash if the comparator is invalid
    function bug(gens::Vector{Int}, exps::Vector{Vector{T}}) where {T}
        inds = collect(1:length(gens))
        cmp  = (x, y) -> Groebner.monom_isless(exps[gens[x]], exps[gens[y]], Groebner._DegRevLex{true}([]))
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
