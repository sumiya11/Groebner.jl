import Primes

@testset "rational reconstruction" begin
    moduli =
        BigInt[107, 199, 509, 2^31 - 1, nextprime(BigInt(2)^100), nextprime(BigInt(2)^200)]
    numbers = [QQ(0), QQ(1), QQ(1, 2), QQ(3, 2), QQ(4, 5), QQ(1, 7), QQ(6)]

    for m in moduli
        bnd = Groebner.ratrec_reconstruction_bound(m)
        for a in numbers
            ac = mod(numerator(a) * invmod(denominator(a), m), m)
            success, (num, den) = Groebner.ratrec_nemo(
                Groebner.Nemo.ZZRingElem(ac),
                Groebner.Nemo.ZZRingElem(m)
            )
            @test Base.unsafe_rational(num, den) == a
        end
    end

    # random checks
    moduli = BigInt[
        17,
        prod(Primes.nextprimes(BigInt(2^10), 20)),
        2^31 - 1,
        prod(Primes.nextprimes(BigInt(2^30 + 3), 5))
    ]
    samples = 1000
    for m in moduli
        bnd = Groebner.ratrec_reconstruction_bound(m)
        for i in 1:samples
            a = BigInt(rand(0:(m - 1)))
            ac = mod(numerator(a) * invmod(denominator(a), m), m)
            success, (num, den) = Groebner.ratrec_nemo(
                Groebner.Nemo.ZZRingElem(ac),
                Groebner.Nemo.ZZRingElem(m)
            )
            if a < sqrt(div(m, 2))
                @test success
            end
            if success
                @test mod(num, m) == mod(a * den, m)
            end
        end
    end
end
