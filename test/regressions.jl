
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

    exps =
        Groebner.ExponentVector{UInt64}[[0x0000, 0x0000], [0x0002, 0x0002], [0x0002, 0x0002]]

    bug(gens, exps)
end
