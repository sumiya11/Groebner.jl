
@testset "packed exponent vector" begin
    PV{T, U} = Groebner.PackedVector{T, U} where {T, U}

    @test Groebner.decompose(0x00030201, UInt8) == (0x01, 0x02, 0x03, 0x00)
    @test Groebner.decompose(0x0a00000400030201, UInt8) == (0x01, 0x02, 0x03, 0x00, 0x04, 0x00, 0x00, 0x0a)

    x = [1, 2, 3, 0, 4]
    ev = Groebner.make_ev(PV{UInt64, UInt8}, x)
    @test ev.data == UInt64[0x0a00000400030201]
    @test typeof(Groebner.totaldeg(ev)) === UInt8
    @test Groebner.totaldeg(ev) == UInt8(10)
    @test Groebner.powertype(ev) === UInt8
    @test Groebner.make_dense(ev, 5) == [1, 2, 3, 0, 4]

    ev = Groebner.make_ev(PV{UInt64, UInt16}, x)
    @test ev.data == UInt64[0x000a000000000004, 0x0000000300020001]
    @test typeof(Groebner.totaldeg(ev)) === UInt16
    @test Groebner.totaldeg(ev) == UInt16(10)
    @test Groebner.powertype(ev) === UInt16

    for B in (UInt8,UInt16,UInt32)
        for T in (UInt64,UInt32,UInt16)
            if sizeof(T) < sizeof(B)
                continue
            end
            for k in 1:10
                n = rand(1:20)
                x = rand(0:10, n)
                ev = Groebner.make_ev(PV{T, B}, x)
                @test x == Groebner.make_dense(ev, n)
            end
        end
    end
end
