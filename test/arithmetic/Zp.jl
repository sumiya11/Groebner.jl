import Primes

@testset "arithmetic in Zp" begin
    m = UInt64(2^30 + 3)
    a = Groebner.DelayedArithmeticZp(UInt64, UInt64, m)
    @test Groebner.n_spare_bits(a) == 1
    @test Groebner.n_safe_consecutive_additions(a) == 2
    @test Groebner.mod_p(UInt64(2^30 + 6), a) === UInt64(3)

    m = UInt64(2^27 + 29)
    a = Groebner.DelayedArithmeticZp(UInt64, UInt32, m)
    @test Groebner.n_spare_bits(a) == 4
    @test Groebner.n_safe_consecutive_additions(a) == 2^7

    m = UInt16(61)
    a = Groebner.DelayedArithmeticZp(UInt16, UInt16, m)
    @test Groebner.n_spare_bits(a) == 2
    @test Groebner.n_safe_consecutive_additions(a) == 8

    for (A, C) in [(UInt64, UInt32), (UInt32, UInt16), (UInt16, UInt8)]
        modulo = Primes.prevprime(typemax(C) >> 2)
        for arithmetic in [
            Groebner.ArithmeticZp(A, C, modulo),
            Groebner.DelayedArithmeticZp(A, C, modulo),
            Groebner.SignedArithmeticZp(signed(A), signed(C), signed(modulo)),
            Groebner.SpecializedArithmeticZp(A, C, modulo)
        ]
            for _ in 1:1000
                T = typeof(Groebner.divisor(arithmetic))
                x = abs(rand(T))
                @test Groebner.mod_p(T(x), arithmetic) == Base.mod(x, T(modulo))
                if gcd(modulo, T(x)) == 1
                    @test Groebner.inv_mod_p(T(x), arithmetic) ==
                          Base.invmod(T(x), T(modulo))
                end
            end
        end
    end
end
