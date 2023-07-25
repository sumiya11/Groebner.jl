import Nemo
import Primes

@testset "Nemo.jl, univariate" begin
    R, x = PolynomialRing(GF(2^31 - 1), "x")
    @test Groebner.groebner([x^2 - 4, x + 2]) == [x + 2]

    R, x = PolynomialRing(QQ, "x")
    @test Groebner.groebner([x^2 - 4, x + 2]) == [x + 2]
end

@testset "Nemo.jl, input-output" begin
    R, (x, y) = AbstractAlgebra.GF(Primes.nextprime(BigInt(2)^100))["x", "y"]
    @test_throws DomainError Groebner.groebner([x, y])

    nemo_orderings_to_test = [:lex, :deglex, :degrevlex]
    nemo_grounds_to_test = [Nemo.GF(2^62 + 135), Nemo.GF(2^31 - 1), Nemo.GF(17), Nemo.QQ]

    for ord in nemo_orderings_to_test
        for ground in nemo_grounds_to_test
            R, x = Nemo.PolynomialRing(ground, "x")
            gb = Groebner.groebner([(x - 1) * (x + 8), (x + 8) * (x + 10)])
            @test gb == [(x + 8)]

            R, (x, y) = Nemo.PolynomialRing(ground, ["x", "y"], ordering=ord)
            fs = [x^2 * y + 3, (2^31 - 5) * x - (2^31 - 4) * y]
            gb = Groebner.groebner(fs)
            @test parent(gb[1]) == R
            @test Groebner.isgroebner(gb)
        end
    end
end
