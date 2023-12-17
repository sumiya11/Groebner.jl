
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

    A, C = UInt64, UInt32
    modulo = C(2^31 - 1)
    for arithmetic in [
        Groebner.ArithmeticZp(A, C, modulo),
        Groebner.DelayedArithmeticZp(A, C, modulo),
        Groebner.SignedArithmeticZp(signed(A), signed(C), signed(modulo)),
        Groebner.SpecializedArithmeticZp(A, C, modulo)
    ]
        for _ in 1:10
            T = typeof(Groebner.divisor(arithmetic))
            x = rand(T)
            @test Groebner.mod_p(T(x), arithmetic) == Base.mod(x, T(modulo))
        end
    end
end
