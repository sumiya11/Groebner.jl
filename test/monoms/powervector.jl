
@testset "power exponent vector" begin
    PV{T} = Groebner.PowerVector{T} where {T}

    x = [1, 2, 3, 0, 4]
    ev = Groebner.make_ev(PV{UInt64}, x)
    @test ev == UInt64.([10, 1, 2, 3, 0, 4])
    @test typeof(Groebner.totaldeg(ev)) === UInt64
    @test Groebner.totaldeg(ev) == UInt64(10)
    @test Groebner.powertype(ev) === Groebner.MonomHash
    tmp = similar(x)
    @test Groebner.make_dense!(tmp, ev) == x
    @test tmp == x

    ev = Groebner.make_ev(PV{UInt32}, x)
    @test ev == UInt32.([10, 1, 2, 3, 0, 4])
    @test typeof(Groebner.totaldeg(ev)) === UInt32
    @test Groebner.totaldeg(ev) == UInt32(10)
    @test Groebner.powertype(ev) === Groebner.MonomHash

    for T in (UInt64, UInt32, UInt16, UInt8)
        for k in 1:10
            n = rand(1:20)
            x = rand(0:10, n)
            ev = Groebner.make_ev(PV{T}, x)
            tmp = similar(x)
            @test x == Groebner.make_dense!(tmp, ev)
            @test x == tmp
            @test Groebner.powertype(ev) === Groebner.MonomHash
        end
    end
end
