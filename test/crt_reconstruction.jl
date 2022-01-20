
using Primes

@testset "Chineese remainder theorem" begin

    modular_images(a, ms) = map(m -> mod(a, m), ms)

    Ms = [
        [2, 3],
        [3, 7],
        [2^31 - 1, Primes.nextprime(2^31)],
        [Primes.nextprime(BigInt(2)^100), Primes.nextprime(BigInt(2)^101)]
    ]
    for ms in Ms
        P = prod(ms)
        nums = [rand(1:P-1) for _ in 1:10]

        for a in nums
            rs = modular_images(a, ms)

            r1, r2 = BigInt.(rs)
            m1, m2 = BigInt.(ms)
            acrt = Groebner.CRT(r1, m1, r2, m2)
            @test acrt == a
        end

    end

end
