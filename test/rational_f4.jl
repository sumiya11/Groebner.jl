
using AbstractAlgebra

@testset "Groebner modular Integers" begin
    R, (x,) = PolynomialRing(ZZ, ["x"], ordering=:degrevlex)
    @test Groebner.groebner([x]) == [x]
    @test Groebner.groebner([4x]) == [x]
    @test Groebner.groebner([x + 1]) == [x + 1]
    @test Groebner.groebner([4x + 1]) == [4x + 1]
    @test Groebner.groebner([x + 5]) == [x + 5]
    @test Groebner.groebner([x + 2^10]) == [x + 2^10]
    @test Groebner.groebner([x + 2^30]) == [x + 2^30]
    @test Groebner.groebner([x + 2^30 + 3]) == [x + 2^30 + 3]
    @test Groebner.groebner([x + 2^31 - 1]) == [x + 2^31 - 1]
    @test Groebner.groebner([x + 2^31 - 1, x^2]) == [1]
    @test Groebner.groebner([3232323x + 7777]) == [3232323x + 7777]
    @test Groebner.groebner([(2^31-1)x + 1]) == [2147483647x + 1]
    @test Groebner.groebner([x + 2147483647]) == [x + 2147483647]

    # consequtive primes: 2^31-1, 2147483659, 2147483693
    # consequtive primes: 2^30+3, 1073741831, 1073741833
    R, (x, y) = PolynomialRing(ZZ, ["x", "y"], ordering=:degrevlex)
    @test Groebner.groebner([4x + 1, 3y + 6]) == [y + 2, 4x + 1]
    @test Groebner.groebner([3y + 6, 4x + 1]) == [y + 2, 4x + 1]
    @test Groebner.groebner([10x + 10y, 10y + 1]) == [10y + 1, 10x - 1]
    @test Groebner.groebner([(2^31-1)x*y + 2147483659x + 2147483693y]) == [(2^31-1)x*y + 2147483659x + 2147483693y]
    @test Groebner.groebner([ (2^30+3)x*y^2 + 1073741831x^2 + 1073741833y + 2^31-1]) == [(2^30+3)x*y^2 + 1073741831x^2 + 1073741833y + 2^31-1]

    Groebner.groebner([4x + 1, 3y + 6])
end

@testset "Groebner modular Rationals" begin
    R, (x,) = PolynomialRing(QQ, ["x"], ordering=:degrevlex)
    @test Groebner.groebner([x]) == [x]
    @test Groebner.groebner([4x]) == [x]
    @test Groebner.groebner([x + 1]) == [x + 1]
    @test Groebner.groebner([4x + 1]) == [x + 1//4]
    @test Groebner.groebner([x + 5]) == [x + 5]
    @test Groebner.groebner([x + 2^10]) == [x + 2^10]
    @test Groebner.groebner([x + 2^30]) == [x + 2^30]
    @test Groebner.groebner([x + 2^30 + 3]) == [x + 2^30 + 3]
    @test Groebner.groebner([x + 2^31 - 1]) == [x + 2^31 - 1]
    @test Groebner.groebner([x + 2^31 - 1, x^2]) == [1]
    @test Groebner.groebner([(3232323//7777)x + 7777//3232323]) == [x + 60481729//10447911976329]
    @test Groebner.groebner([((2^31-1)//1)x + 1]) == [x + 1//2147483647]
    @test Groebner.groebner([(1//(2^31-1))x + 1]) == [x + 2147483647]
    @test Groebner.groebner([1//(2^30+3)*x^2 + (2^30+3)x + 1//(1073741831)]) == [x^2 + 1152921511049297929x + (2^30+3)//1073741831]

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)

    fs = [x, y]
    G = Groebner.groebner(fs)
    @test G ≂ [y, x]

    fs = [5y, 3//4*x]
    G = Groebner.groebner(fs)
    @test G ≂ [y, x]

    fs = [x^2 - (1//6)*y, x*y]
    G = Groebner.groebner(fs)
    @test G ≂ [y^2, x*y, x^2 - (1//6)y]

    fs = [QQ(11, 3)*x^2*y - QQ(2, 4)*x*y - QQ(1, 7)*y,
            QQ(1, 1)*x*y + QQ(7, 13)*x]
    G = Groebner.groebner(fs)
    @test G ≂ [y^2 + 7//13*y,
                x*y + 7//13*x,
                x^2 - 3//22*x + 39//539*y]

    
    root = Groebner.change_ordering(Groebner.rootn(3, ground=QQ), :degrevlex)
    gb = Groebner.groebner(root)
    @test Groebner.isgroebner(gb)

    root = Groebner.change_ordering(Groebner.rootn(4, ground=QQ), :degrevlex)
    gb = Groebner.groebner(root)
    @test Groebner.isgroebner(gb)

    root = Groebner.change_ordering(Groebner.rootn(8, ground=QQ), :degrevlex)
    gb = Groebner.groebner(root)
    @test Groebner.isgroebner(gb)
    

    
    noon = Groebner.change_ordering(Groebner.noonn(3, ground=QQ), :degrevlex)
    gb = Groebner.groebner(noon)
    @test Groebner.isgroebner(gb)

    noon = Groebner.change_ordering(Groebner.noonn(4, ground=QQ), :degrevlex)
    gb = Groebner.groebner(noon)
    @test Groebner.isgroebner(gb)
    
end

@testset "Groebner modular corner cases" begin

    R, (x, y, z, w) = PolynomialRing(QQ, ["x", "y", "z", "w"], ordering=:degrevlex)

    fs = [(12345678//12347)x,
          (222222221111123//2131232232097)y + z]
    G = Groebner.groebner(fs)
    @test G == [y + 2131232232097//222222221111123*z, x]

    fs = [
           (224//225)*x + y,
        x^2 + 1]
    ans = [x + 225//224*y, y^2 + 50176//50625]
    @test Groebner.groebner(fs) == ans


    fs = [
        (12345//12345678)x + y,
        x^2 + 1]
    ans = [x + 4115226//4115*y, y^2 + 16933225//16935085031076]
    @test Groebner.groebner(fs) == ans

    fs = [
        (12345//12345678)x + y,
        x^2 + z^2,
        y^2 + w^2,
        x + 2*y + 3*z + 4*w
    ]
    G = Groebner.groebner(fs)

    @test G ≂ [w^3,
                z*w + 1692834523553//4063974*w^2,
                z^2 - 16935085031076//16933225*w^2,
                y - 12345//4106996*z - 4115//1026749*w,
                x + 6172839//2053498*z + 4115226//1026749*w]

    # what if we are unlucky to start from initial prime in coefficients
    G = [
        (2^30 + 3)x,
        (1//(2^30 + 3))y
    ]
    G = Groebner.groebner(G)
    @test G == [y, x]

    G = [
        (2^31 - 1)*x + y,
        (1//(2^31 - 1))*y + 1
    ]
    G = Groebner.groebner(G)
    @test G == [y + 2147483647, x - 1]

    R, (x, y, z, w) = PolynomialRing(QQ, ["x", "y", "z", "w"], ordering=:lex)
    fs = [(12345678//12347)x,
          (222222221111123//2131232232097)y + z]
    G = Groebner.groebner(fs)
    @test G == [y + 2131232232097//222222221111123*z, x]
end
