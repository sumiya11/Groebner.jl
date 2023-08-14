using Primes

@testset "rational reconstruction" begin
    # these are just small sanity checks 
    moduli =
        BigInt[107, 199, 509, 2^31 - 1, nextprime(BigInt(2)^100), nextprime(BigInt(2)^200)]
    numbers = [QQ(0), QQ(1), QQ(1, 2), QQ(3, 2), QQ(4, 5), QQ(1, 7), QQ(6)]

    num, den, buf, buf1, buf2, buf3, u1, u2, u3, v1, v2, v3 = [BigInt(0) for _ in 1:13]
    for m in moduli
        bnd = Groebner.rational_reconstruction_bound(m)
        for a in numbers
            ac = numerator(a) * invmod(denominator(a), m)
            ar1 = Groebner.rational_reconstruction!(
                num,
                den,
                bnd,
                buf,
                buf1,
                buf2,
                buf3,
                u1,
                u2,
                u3,
                v1,
                v2,
                v3,
                ac,
                m
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
        bnd = Groebner.rational_reconstruction_bound(m)
        for i in 1:samples
            a = BigInt(rand(0:(m - 1)))
            success = Groebner.rational_reconstruction!(
                num,
                den,
                bnd,
                buf,
                buf1,
                buf2,
                buf3,
                u1,
                u2,
                u3,
                v1,
                v2,
                v3,
                a,
                m
            )
            if a < sqrt(m)
                @test success
            end
            if success
                @test mod(num, m) == mod(a * den, m)
            end
        end
    end
end
