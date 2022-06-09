
@testset "fglm finite field" begin

    R, (x, y) = PolynomialRing(GF(2^31-1), ["x","y"], ordering=:degrevlex)

    noon = Groebner.change_ordering(Groebner.noonn(2, ground=GF(2^31-1)), :degrevlex)
    gb = Groebner.groebner(noon, ordering=:degrevlex)
    x1, x2 = gens(parent(first(Groebner.fglm(gb))))
    @test Groebner.fglm(gb) == [x2^5 + 1952257860*x2^4 + 1288490186*x2^3 + 2*x2^2 + 1028839893*x2 + 644245093,
                                x1 + x2^4 + 1952257860*x2^3 + 644245093*x2^2 + x2 + 1952257860]

    gb = Groebner.groebner(noon, ordering=:deglex)
    x1, x2 = gens(parent(first(Groebner.fglm(gb))))
    @test Groebner.fglm(gb) == [x2^5 + 1952257860*x2^4 + 1288490186*x2^3 + 2*x2^2 + 1028839893*x2 + 644245093,
                                x1 + x2^4 + 1952257860*x2^3 + 644245093*x2^2 + x2 + 1952257860]

    gb = Groebner.groebner(noon, ordering=:lex)
    x1, x2 = gens(parent(first(Groebner.fglm(gb))))
    @test Groebner.fglm(gb) == [x2^5 + 1952257860*x2^4 + 1288490186*x2^3 + 2*x2^2 + 1028839893*x2 + 644245093,
                                x1 + x2^4 + 1952257860*x2^3 + 644245093*x2^2 + x2 + 1952257860]

end

@testset "fglm rationals" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x","y","z"], ordering=:degrevlex)

    gb = [  z^2 - (9//490)*y - 201//980*z + 13//980,
            y*z - 4//35*y + 2//35*z - 1//35,
            y^2 - 1//10*y - 6//5*z + 1//10,
            x^2 - 4//5*y + 2//5*z - 1//5]

    x, y, z = gens(parent(first(Groebner.fglm(gb))))
    @test Groebner.fglm(gb) == [
                                z^3 - 313//980*z^2 + 37//980*z - 1//490,
                                y - 490//9*z^2 + 67//6*z - 13//18,
                                x^2 - 392//9*z^2 + 28//3*z - 7//9,
                                ]

end
