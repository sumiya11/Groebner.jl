using Primes

@testset "chinese remainder theorem" begin
    modular_images(a, ms) = map(m -> mod(a, m), ms)

    Ms = [
        [2, 3],
        [3, 7],
        [2^31 - 1, Primes.nextprime(2^31)],
        [Primes.nextprime(BigInt(2)^100), Primes.nextprime(BigInt(2)^31)],
        [Primes.nextprime(BigInt(2)^2000), Primes.nextprime(BigInt(2)^31)]
    ]
    buf, buf1, buf2, buf3, n1, n2 = [BigInt(0) for _ in 1:6]

    for (m1, m2) in Ms
        M = prod((m1, m2))
        nums = [rand(0:(M - 1)) for _ in 1:100]
        minv1, minv2 = invmod(m1, m2), invmod(m2, m1)
        for a in nums
            as = modular_images(a, (m1, m2))
            a1, a2 = BigInt.(as)
            m1, m2 = BigInt.((m1, m2))
            minv1, minv2 = BigInt.((minv1, minv2))
            c1, c2 = m2 * minv2, m1 * minv1
            Groebner.CRT!(M, buf, n1, n2, a1, c1, UInt64(a2), c2)
            @test buf == a
        end
    end
end
