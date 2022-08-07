
using AbstractAlgebra

@testset "groebner zeros" begin
    for field in [GF(2), GF(2^31-1), ZZ, QQ]
        R, (x, y, z) = PolynomialRing(field, ["x","y","z"], ordering=:lex)

        @test Groebner.groebner([x - x]) == [R(0)]
        @test Groebner.groebner([R(0), R(0), R(0)]) == Groebner.groebner([R(0)]) == [R(0)]
        @test Groebner.groebner([x, R(0), y, R(0)]) == [y, x]
        
        @test_throws DomainError Groebner.groebner([])
    end
end

@testset "isgroebner zeros" begin
    for field in [GF(2), GF(2^31-1), ZZ, QQ]
        R, (x, y, z) = PolynomialRing(field, ["x","y","z"], ordering=:lex)

        @test Groebner.isgroebner([x - x])
        @test Groebner.isgroebner([R(0), R(0), R(0)])
        @test Groebner.isgroebner([x, R(0), y, R(0)])
        
        @test_throws DomainError Groebner.isgroebner([])
    end
end

@testset "normalform zeros" begin
    for field in [GF(2), GF(2^31-1), ZZ, QQ]
        R, (x, y, z) = PolynomialRing(field, ["x","y","z"], ordering=:lex)

        @test Groebner.normalform([R(0)], [R(0)]) == [R(0)]
        @test Groebner.normalform([R(0), R(0), R(0)], x) == x
        # @test Groebner.normalform([R(0), R(0), R(0)], R(0)) == R(0)
        
        # @test_throws DomainError Groebner.normalform([], x)
    end
end

@testset "fglm zeros" begin
    for field in [GF(2), GF(2^31-1), ZZ, QQ]
        R, (x, y, z) = PolynomialRing(field, ["x","y","z"], ordering=:degrevlex)

        b = Groebner.fglm([R(0)])
        @test b == [zero(parent(first(b)))]
        b = Groebner.fglm([R(0), R(0), R(0)])
        @test b == [zero(parent(first(b)))]
        
        @test_throws DomainError Groebner.fglm([])
    end
end

@testset "kbase zeros" begin
    for field in [GF(2), GF(2^31-1), ZZ, QQ]
        R, (x, y, z) = PolynomialRing(field, ["x","y","z"], ordering=:lex)

        @test_throws DomainError Groebner.kbase([R(0)]) 
        @test_throws DomainError Groebner.kbase([R(0), R(0), R(0)])
    
        @test Groebner.kbase([x, y^2, R(0), z, R(0)]) == [R(1), y]
    end
end