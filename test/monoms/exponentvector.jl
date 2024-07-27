
@testset "exponent vector" begin
    PV{T} = Groebner.ExponentVector{T} where {T}

    x = [1, 2, 3, 0, 4]
    ev = Groebner.monom_construct_from_vector(PV{UInt64}, x)
    @test ev == UInt64.([10, 1, 2, 3, 0, 4])
    @test typeof(Groebner.monom_totaldeg(ev)) === UInt64
    @test Groebner.monom_totaldeg(ev) == UInt64(10)
    @test Groebner.monom_entrytype(ev) === UInt64
    tmp = similar(x)
    @test Groebner.monom_to_vector!(tmp, ev) == x
    @test tmp == x

    @test Groebner._monom_overflow_check(ev)

    ev = Groebner.monom_construct_from_vector(PV{UInt32}, x)
    @test ev == UInt32.([10, 1, 2, 3, 0, 4])
    @test typeof(Groebner.monom_totaldeg(ev)) === UInt32
    @test Groebner.monom_totaldeg(ev) == UInt32(10)
    @test Groebner.monom_entrytype(ev) === UInt32

    for T in (UInt64, UInt32, UInt16, UInt8)
        for k in 1:10
            n = rand(1:20)
            x = rand(0:6, n)
            ev = Groebner.monom_construct_from_vector(PV{T}, x)
            tmp = similar(x)
            @test x == Groebner.monom_to_vector!(tmp, ev)
            @test x == tmp
            @test Groebner.monom_entrytype(ev) === T
        end
    end
end
