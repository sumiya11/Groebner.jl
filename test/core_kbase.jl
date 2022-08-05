
@testset "kbase simple" begin
    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)

    @test Groebner.kbase([x, y]) == [R(1)]
    @test Groebner.kbase([x + y, y]) == [R(1)]
    @test Groebner.kbase([-2x + 8, 10y + 1]) == [R(1)]

    @test Groebner.kbase([x^2 + 1, y^2 - x - 1]) == [R(1), y, x, x*y]
    @test Groebner.kbase([8x^5 + 1, y - 1]) == [R(1), x, x^2, x^3, x^4]
    @test Groebner.kbase([x^2, y^2, x*y]) == [R(1), y, x]
    @test Groebner.kbase([x^3, y^3, x*y]) == [R(1), y, x, y^2, x^2]
    @test length(Groebner.kbase([8x^77 + 1, y^34 - 1])) == 77*34


    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:lex)
    @test Groebner.kbase([x, y]) == [R(1)]
    @test Groebner.kbase([x + y, y]) == [R(1)]
    @test Groebner.kbase([-2x + 8, 10y + 1]) == [R(1)]

    @test Groebner.kbase([x^2 + 1, y^2 - 1]) == [R(1), y, x, x*y]
    @test Groebner.kbase([8x^5 + 1, y - 1]) == [R(1), x, x^2, x^3, x^4]
    @test Groebner.kbase([x^2, y^2, x*y]) == [R(1), y, x]
    @test Groebner.kbase([x^3, y^3, x*y]) == [R(1), y, y^2, x, x^2]
    @test length(Groebner.kbase([8x^77 + 1, y^34 - 1])) == 77*34


    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:deglex)

    @test Groebner.kbase([x, y]) == [R(1)]
    @test Groebner.kbase([x + y, y]) == [R(1)]
    @test Groebner.kbase([-2x + 8, 10y + 1]) == [R(1)]

    @test Groebner.kbase([x^2 + 1, y^2 - x - 1]) == [R(1), y, x, x*y]
    @test Groebner.kbase([8x^5 + 1, y - 1]) == [R(1), x, x^2, x^3, x^4]
    @test Groebner.kbase([x^2, y^2, x*y]) == [R(1), y, x]
    @test Groebner.kbase([x^3, y^3, x*y]) == [R(1), y, x, y^2, x^2]
    @test length(Groebner.kbase([8x^77 + 1, y^34 - 1])) == 77*34
end

@testset "kbase big" begin
    root4 = Groebner.rootn(4)
    gb = Groebner.groebner(root4)
    R = parent(first(root4))
    (x1,x2,x3,x4) = gens(R)
    @test Groebner.kbase(gb) == [R(1), x4, x4^2, x4^3, x3, x3*x4, x3*x4^2,
                                 x3*x4^3, x3^2, x3^2*x4, x3^2*x4^2, x3^2*x4^3,
                                 x2, x2*x4, x2*x4^2, x2*x4^3, x2*x3, x2*x3*x4,
                                 x2*x3*x4^2, x2*x3*x4^3, x2*x3^2, x2*x3^2*x4,
                                 x2*x3^2*x4^2, x2*x3^2*x4^3]

     root6 = Groebner.rootn(6)
     gb = Groebner.groebner(root6)
     @test length(Groebner.kbase(gb)) == 720
     root6 = Groebner.rootn(6)
     gb = Groebner.groebner(root6, ordering=:deglex)
     @test length(Groebner.kbase(gb)) == 720
     root6 = Groebner.rootn(6)
     gb = Groebner.groebner(root6, ordering=:degrevlex)
     @test length(Groebner.kbase(gb)) == 720

     noon2 = Groebner.change_ordering(Groebner.noonn(2), :degrevlex)
     R = parent(first(noon2))
     (x1,x2) = gens(R)
     gb = Groebner.groebner(noon2, ordering=:degrevlex)
     @test Groebner.kbase(gb) == [R(1), x2, x1, x2^2, x1*x2]

     noon3 = Groebner.noonn(3)
     R = parent(first(noon3))
     gb = Groebner.groebner(noon3, ordering=:degrevlex)
     @test length(Groebner.kbase(gb)) == 21

     noon3 = Groebner.noonn(7)
     R = parent(first(noon3))
     gb = Groebner.groebner(noon3, ordering=:degrevlex)
     @test length(Groebner.kbase(gb)) == 2173
end
