
@testset "groebner orderings" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x","y","z"], ordering=:deglex)

    Groebner.groebner([x, y], ordering=Groebner.Lex())
    Groebner.groebner([x, y], ordering=Groebner.DegLex())
    Groebner.groebner([x, y], ordering=Groebner.DegRevLex())
    
    Groebner.groebner([x, y], ordering=Groebner.WeightedOrdering([1, 2, 3]))
    @test_throws DomainError Groebner.groebner([x, y], ordering=Groebner.WeightedOrdering([1, 0]))
    @test_throws DomainError Groebner.groebner([x, y], ordering=Groebner.WeightedOrdering([1, 0, 1, 9]))
    
    ord = Groebner.BlockOrdering(1:1, Groebner.WeightedOrdering([0]), 2:3, Groebner.WeightedOrdering([2, 3]))
    Groebner.groebner([x, y], ordering=ord)
    ord = Groebner.BlockOrdering(1:1, Groebner.WeightedOrdering([0, 0]), 2:3, Groebner.WeightedOrdering([2, 3]))
    @test_throws DomainError Groebner.groebner([x, y], ordering=ord)
    ord = Groebner.BlockOrdering(1:1, Groebner.WeightedOrdering([0]), 2:3, Groebner.WeightedOrdering([2]))
    @test_throws DomainError Groebner.groebner([x, y], ordering=ord)

    R, (x1,x2,x3,x4,x5,x6) = QQ["x1","x2","x3","x4","x5","x6"];
    ord_3_6 = Groebner.BlockOrdering(3:4, Groebner.DegRevLex(), 5:6, Groebner.WeightedOrdering([0, 3]))
    ord = Groebner.BlockOrdering(1:2, Groebner.Lex(), 3:6, ord_3_6)
    Groebner.groebner([x1, x2], ordering=ord)

    ord_3_6 = Groebner.BlockOrdering(3:3, Groebner.WeightedOrdering([0]), 4:6, Groebner.WeightedOrdering([1, 0, 3]))
    ord = Groebner.BlockOrdering(1:2, Groebner.WeightedOrdering([1, 1]), 3:6, ord_3_6)
    Groebner.groebner([x1, x2], ordering=ord)

    ord_3_6 = Groebner.BlockOrdering(3:4, Groebner.WeightedOrdering([0, 1, 2]), 5:6, Groebner.WeightedOrdering([1, 0]))
    ord = Groebner.BlockOrdering(1:2, Groebner.WeightedOrdering([1, 1]), 3:6, ord_3_6)
    @test_throws DomainError Groebner.groebner([x1, x2], ordering=ord)

    ord = Groebner.MatrixOrdering([
        1 0 0 0  1  2;
        0 1 0 0 -2 -1;
        0 0 1 0  0  0;
        0 0 0 1  0  0;
    ])
    Groebner.groebner([x1, x2], ordering=ord)

    ord = Groebner.MatrixOrdering([
        1 0 0 0;
        0 1 0 0;
        0 0 1 0;
    ])
    @test_throws DomainError Groebner.groebner([x1, x2], ordering=ord)
end
