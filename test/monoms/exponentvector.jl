
@testset "exponent vector" begin
    PV{T} = Groebner.ExponentVector{T} where {T}

    x = [1, 2, 3, 0, 4]
    ev = Groebner.construct_monom(PV{UInt64}, x)
    @test ev == UInt64.([10, 1, 2, 3, 0, 4])
    @test typeof(Groebner.totaldeg(ev)) === UInt64
    @test Groebner.totaldeg(ev) == UInt64(10)
    @test Groebner.entrytype(ev) === Groebner.MonomHash
    tmp = similar(x)
    @test Groebner.monom_to_dense_vector!(tmp, ev) == x
    @test tmp == x

    @test Groebner._monom_overflow_check(ev)

    ev = Groebner.construct_monom(PV{UInt32}, x)
    @test ev == UInt32.([10, 1, 2, 3, 0, 4])
    @test typeof(Groebner.totaldeg(ev)) === UInt32
    @test Groebner.totaldeg(ev) == UInt32(10)
    @test Groebner.entrytype(ev) === Groebner.MonomHash

    for T in (UInt64, UInt32, UInt16, UInt8)
        for k in 1:10
            n = rand(1:20)
            x = rand(0:10, n)
            ev = Groebner.construct_monom(PV{T}, x)
            tmp = similar(x)
            @test x == Groebner.monom_to_dense_vector!(tmp, ev)
            @test x == tmp
            @test Groebner.entrytype(ev) === Groebner.MonomHash
        end
    end
end
