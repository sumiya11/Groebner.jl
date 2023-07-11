using Primes

@testset "chinese remainder theorem" begin
    modular_images(a, ms) = map(m -> mod(a, m), ms)

    Ms = [
        [2, 3],
        [3, 7],
        [2^31 - 1, Primes.nextprime(2^31)],
        [Primes.nextprime(BigInt(2)^100), Primes.nextprime(BigInt(2)^31)],
        [Primes.nextprime(BigInt(2)^10000), Primes.nextprime(BigInt(2)^31)]
    ]
    buf, buf1, buf2, buf3, n1, n2 = [BigInt(0) for _ in 1:6]

    for ms in Ms
        M = prod(ms)
        nums = [rand(0:(M - 1)) for _ in 1:1000]

        msinv = [gcdx(ms[1], ms[2])[2], gcdx(ms[1], ms[2])[3]]
        for a in nums
            rs = modular_images(a, ms)

            r1, r2 = BigInt.(rs)
            m1, m2 = BigInt.(ms)
            minv1, minv2 = BigInt.(msinv)
            Groebner.CRT!(M, buf, n1, n2, r1, minv1, UInt64(r2), minv2, m1, m2)
            @test buf == a
        end
    end
end
