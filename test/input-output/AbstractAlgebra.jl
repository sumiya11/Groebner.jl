using AbstractAlgebra

@testset "AbstractAlgebra.jl, univariate" begin
    R, x = PolynomialRing(GF(2^31 - 1), "x")
    @test Groebner.groebner([x^2 - 4, x + 2]) == [x + 2]

    R, x = PolynomialRing(QQ, "x")
    @test Groebner.groebner([x^2 - 4, x + 2]) == [x + 2]

    R, x = PolynomialRing(ZZ, "x")
    @test Groebner.groebner([x^2 - 4, x + 2]) == [x + 2]
end

@testset "AbstractAlgebra.jl, input-output" begin
    # TODO: gracefully error if the input is not supported
    # e.g., 
    # R, (x, y) = AbstractAlgebra.GF(nextprime(BigInt(2)^100))["x","y"]
    
    aa_orderings_to_test = [:lex, :degrevlex, :deglex]
    aa_grounds_to_test = [AbstractAlgebra.GF(2^31 - 1), AbstractAlgebra.QQ]

    for ord in aa_orderings_to_test
        for ground in aa_grounds_to_test
            R, (x,) = PolynomialRing(ground, ["x"], ordering=ord)
            @test parent(first(Groebner.groebner([x]))) == R

            R, (x, y) = PolynomialRing(ground, ["x", "y"], ordering=ord)
            fs = [x^2 * y + 3, (2^31 - 5) * x - (2^31 - 4) * y]
            gb = Groebner.groebner(fs)
            @test parent(gb[1]) == R
            @test all(
                in([
                    y^3 + ground(13835057990857654347) // ground(4611686001247518736),
                    x - ground(2147483644) // ground(2147483643) * y
                ]),
                gb
            )

            @test Groebner.isgroebner(gb)
            @test all(iszero, Groebner.normalform(gb, gb))
        end
    end
end
