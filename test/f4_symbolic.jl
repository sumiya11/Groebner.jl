
using .GroebnerBases: selectnormal, leftright, symbolic_preprocessing,
                    pairdegree, update!, reduction

R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])


@testset "F4 pairdegree" begin
    @test pairdegree(x, x) == 1
    @test pairdegree(x, R(1)) == 1
    @test pairdegree(x, y*z) == 3
    @test pairdegree(x^2*y + 1, y*z^2 + z^5) == 5
end

@testset "F4 selection" begin

    fs = [ (x*y, y*z) ]
    selected, not = copy(fs), similar(fs, 0)
    @test selectnormal(fs) == selected

    fs = [ (x, R(1)), (R(1), z) ]
    selected, not = copy(fs), similar(fs, 0)
    @test selectnormal(fs) == selected

    fs = [ (x, R(1)), (R(1), x) ]
    selected = selectnormal(copy(fs))
    @test length(selected) == 2 && fs == selected

    fs = [ (R(5), 10x), (R(3), 2y), (x*y, x) ]
    selected = selectnormal(copy(fs))
    @test selected ≂ fs[1:2] # && not == [fs[end]]

    fs = [
        (x*y*z + x, z^2), (x^4, R(8)), (x^2 + z, y^2 - y), (1 - z^2, x^2 + x*y),
        (R(1), x^5), (x^4, y^4), (x^2 + 2x*y + y^2, z^3 + z + 1),
        (10x*z - R(5), (1//100)z^2*y)
    ]
    selected = union(fs[1:4], [fs[end]])
    not = filter(f -> !(f in selected), fs)
    @test selectnormal(copy(fs)) ≂ selected

end

@testset "F4 pair construction" begin

    # the leftright function is not very typestable yet
    T = typeof((x*y, R(-1)))
    fs = [ (x*y, R(-1)) ]
    @test leftright(fs) ≂ T[ (R(1), x*y), (x*y, R(-1)) ]

    fs = [ (x*y, x), (x, x*y) ]
    @test leftright(fs) ≂ T[ (y, x), (R(1), x*y) ]

    fs = [ (5x + 1, -4y + 3z) ]
    @test leftright(fs) ≂ T[ (y, 5x + 1), (x, -4y + 3z) ]

end

@testset "F4 symbolic preprocessing" begin

    G = [x*y^2, y, x + 1, x^2]
    Ld = [(x, y), (x, x + 1), (y, x + 1), (1, x^2)]
    @test symbolic_preprocessing(Ld, G, []) ⊂ [
        x*y, x^2 + x, x*y + y, x^2, x + 1, y
    ]

    G = [R(1), x, x^2, x^3]
    Ld = [(x, 1), (1, x)]
    @test symbolic_preprocessing(Ld, G, []) ⊂ [x]

    G = [x, y, z]
    Ld = [(y, x), (z, x), (z, y), (x, y), (x, z), (y, z)]
    # We need to decide on the the uniqueness of the polynomial
    @test symbolic_preprocessing(Ld, G, []) ⊂ [x*y, x*z, y*z]

    G = [y^2, x*y + y^4, R(1)]
    Ld = [(1, y^2), (1, x*y + y^4), (y^2, 1), (x*y, 1)]
    @test symbolic_preprocessing(Ld, G, []) ⊂ [
        y^2, x*y + y^4, x*y, y^4
    ]

    G = [5y^2, x*y + y^4, R(1)]
    Ld = [(1, 5*y^2), (1, x*y + y^4), (y^2, 1), (x*y, 1)]
    # We need to decide on the the uniqueness of the polynomial
    @test symbolic_preprocessing(Ld, G, []) ⊂ [
        5*y^2, x*y + y^4, x*y, 5*y^4, y^4, y^2
    ]
end

@testset "F4 reduction" begin


end

@testset "F4 update" begin


end
