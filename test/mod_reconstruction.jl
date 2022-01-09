
using Primes


@testset "Rational reconstruction" begin
    ms = [107, 199, 509, 2^31 - 1]
    as = [QQ(0), QQ(1), QQ(1, 2), QQ(3, 2), QQ(4, 5), QQ(1, 7), QQ(6)]

    for m in ms
        for a in as
            gf = GF(m)
            ac = lift( numerator(a) // gf(denominator(a)) )
            ar1 = GroebnerBases.rational_reconstruction(ac, BigInt(m))
            ar2 = GroebnerBases.rational_reconstruction(Int(ac), m)
            @test ar1 == ar2 == a
        end
    end

    # random tests
    ms = Primes.primes(2^20, 2^21)[1:3]
    as = [QQ(rand(1:100), rand(1:100)) for _ in 1:10]
    for m in ms
        for a in as
            gf = GF(m)
            ac = lift( numerator(a) // gf(denominator(a)) )
            ar1 = GroebnerBases.rational_reconstruction(ac, BigInt(m))
            ar2 = GroebnerBases.rational_reconstruction(Int(ac), m)
            @test ar1 == ar2 == a
        end
    end

end
