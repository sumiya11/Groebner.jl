
@testset "arithmetic in Zp" begin
    m = UInt64(2^30 + 3)
    a = Groebner.UseSpareBitsArithmeticZp(m)
    @test Groebner.n_spare_bits(a) == 1
    @test Groebner.n_safe_consecutive_additions(a) == 2
    @test Groebner.mod_p(UInt64(2^30 + 6), a) === UInt64(3)

    m = UInt64(2^27 + 29)
    a = Groebner.UseSpareBitsArithmeticZp(m)
    @test Groebner.n_spare_bits(a) == 4
    @test Groebner.n_safe_consecutive_additions(a) == 2^7

    m = UInt16(61)
    a = Groebner.UseSpareBitsArithmeticZp(m)
    @test Groebner.n_spare_bits(a) == 2
    @test Groebner.n_safe_consecutive_additions(a) == 8

    m = UInt64(2^31 + 11)
    a = Groebner.UseSpareBitsArithmeticZp(m)
    @test Groebner.n_spare_bits(a) == 0
    @test Groebner.n_safe_consecutive_additions(a) == 0
end
