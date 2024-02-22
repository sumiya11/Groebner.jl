import Primes

@testset "arithmetic in Zp" begin
    a = Groebner.select_arithmetic(UInt32, 2^31 - 1, :auto, false)
    @test a isa Groebner.SpecializedArithmeticZp{UInt64, UInt32}
    a = Groebner.select_arithmetic(UInt64, 2^31 - 1, :auto, true)
    @test a isa Groebner.SpecializedArithmeticZp{UInt64, UInt64}

    a = Groebner.select_arithmetic(UInt32, Primes.prevprime(2^27 - 1), :auto, false)
    @test a isa Groebner.DelayedArithmeticZp{UInt64, UInt32}
    a = Groebner.select_arithmetic(UInt64, Primes.prevprime(2^27 - 1), :auto, true)
    @test a isa Groebner.DelayedArithmeticZp{UInt64, UInt64}
    a = Groebner.select_arithmetic(UInt64, Primes.prevprime(2^28 - 1), :auto, true)
    @test a isa Groebner.SpecializedArithmeticZp{UInt64, UInt64}
    a = Groebner.select_arithmetic(UInt32, Primes.prevprime(2^28 - 1), :auto, false)
    @test a isa Groebner.SpecializedArithmeticZp{UInt64, UInt32}

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

@testset "arithmetic in Zp x 4" begin
    m = (Primes.nextprimes(2^30 + 3, 4)...,)
    i32_4x = Groebner.CompositeInt{4, Int32}
    i64_4x = Groebner.CompositeInt{4, Int64}
    m_4x = i32_4x(m)
    arithmetic = Groebner.SignedCompositeArithmeticZp(i64_4x, i32_4x, m_4x)

    for _ in 1:1000
        T = Int64
        x = ntuple(_ -> rand(T), 4)
        @test Groebner.mod_p(i64_4x(x), arithmetic).data == Base.mod.(x, m)
        if all(==(1), gcd.(m, x))
            @test Groebner.inv_mod_p(i64_4x(x), arithmetic).data == Base.invmod.(x, m)
        end
    end
end
