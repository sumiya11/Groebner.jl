

R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])

using .GroebnerBases: selectnormal, leftright, symbolic_preprocessing,
                    pairdegree, update!

@testset "F4 pairdegree" begin

end

@testset "F4 selection" begin

    fs = [ (x*y, y*z) ]
    selected, not = copy(fs), similar(fs, 0)
    @test selectnormal(fs) == ( selected, not )

    fs = [ (x, R(1)), (R(1), z) ]
    selected, not = copy(fs), similar(fs, 0)
    @test selectnormal(fs) == ( selected, not )

    fs = [ (x, R(1)), (R(1), x) ]
    selected, not = selectnormal(copy(fs))
    @test length(selected) == 2 && all(x -> x in fs, selected) && isempty(not)

    fs = [ (R(5), 10x), (R(3), 2y), (x*y, x) ]
    selected, not = selectnormal(copy(fs))
    @test Set(selected) == Set(fs[1:2]) && not == [fs[end]]

    fs = [
        (x*y*z + x, z^2), (x^4, R(8)), (x^2 + z, y^2 - y), (1 - z^2, x^2 + x*y),
        (R(1), x^5), (x^4, y^4), (x^2 + 2x*y + y^2, z^3 + z + 1),
        (10x*z - R(5), (1//100)z^2*y)
    ]
    selected = union(fs[1:4], [fs[end]])
    not = filter(f -> !(f in selected), fs)
    @test Set.((selectnormal(copy(fs)))) == Set.(( selected, not ))

end

@testset "F4 pair construction" begin

    # the leftright function is not very typestable yet
    T = typeof((x*y, R(-1)))
    fs = [ (x*y, R(-1)) ]
    @test Set{T}( leftright(fs) ) == Set([ (R(1), x*y), (x*y, R(-1)) ])

    fs = [ (x*y, x), (x, x*y) ]
    @test Set{T}( leftright(fs) ) == Set([ (y, x), (R(1), x*y) ])

    fs = [ (5x + 1, -4y + 3z) ]
    @test Set{T}( leftright(fs) ) == Set([ (y, 5x + 1), (x, -4y + 3z) ])

end

@testset "F4 symbolic preprocessing" begin


end

@testset "F4 updating" begin


end
