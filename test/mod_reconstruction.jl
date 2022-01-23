
using Primes


@testset "Rational reconstruction" begin
    
    ms = BigInt[107, 199, 509, 2^31 - 1, nextprime(BigInt(2)^100), nextprime(BigInt(2)^200)]
    as = [QQ(0), QQ(1), QQ(1, 2), QQ(3, 2), QQ(4, 5), QQ(1, 7), QQ(6)]

    num, den, buf, buf1, buf2, buf3, u1, u2, u3, v1, v2, v3 = [BigInt(0) for _ in 1:13]
    for m in ms
        bnd = Groebner.rational_bound(m)
        for a in as
            ac = numerator(a) * invmod(denominator(a), m)
            ar1 = Groebner.rational_reconstruction!(
                        num, den, bnd, buf, buf1, buf2, buf3,
                        u1, u2, u3, v1, v2, v3, ac, m)
            ar2 = Groebner.rational_reconstruction(ac, m)
            @test Base.unsafe_rational(num, den) == ar2 == a
        end
    end

    # random tests
    ms = Primes.nextprimes(BigInt(2^31-1), 100)
    ases = 100
    for m in ms
        bnd = Groebner.rational_bound(m)
        for i in 1:ases
            ac = rand(0:m-1)
            ar1 = Groebner.rational_reconstruction(ac, m)
            ar2 = Groebner.rational_reconstruction!(
                        num, den, bnd, buf, buf1, buf2, buf3,
                        u1, u2, u3, v1, v2, v3, ac, m)
            @test Base.unsafe_rational(num, den) == ar1
        end
    end

    # Big tests
    ms = Primes.nextprimes(BigInt(2)^100, 10)
    ases = 100
    for m in ms
        bnd = Groebner.rational_bound(m)
        for i in 1:ases
            ac = rand(0:m-1)
            ar1 = Groebner.rational_reconstruction(ac, m)
            ar2 = Groebner.rational_reconstruction!(
                        num, den, bnd, buf, buf1, buf2, buf3,
                        u1, u2, u3, v1, v2, v3, ac, m)
            @test Base.unsafe_rational(num, den) == ar1
        end
    end

end
