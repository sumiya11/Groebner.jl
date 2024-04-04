
@testset "sparse exponent vector" begin
    PV{T, I, N} = Groebner.SparseExponentVector{T, I, N} where {T, I, N}

    ev = Groebner.monom_construct_const(PV{UInt16, Int, 60}, 60)
    @test Groebner.monom_totaldeg(ev) == UInt16(0)
    @test Groebner.monom_entrytype(ev) === Groebner.MonomHash

    x = [1, 2, 3, 0, 4]
    ev = Groebner.monom_construct_from_vector(PV{UInt64, Int, 5}, x)
    @test typeof(Groebner.monom_totaldeg(ev)) === UInt64
    @test Groebner.monom_totaldeg(ev) == UInt64(10)
    @test Groebner.monom_entrytype(ev) === Groebner.MonomHash
    tmp = similar(x)
    @test Groebner.monom_to_vector!(tmp, ev) == x
    @test tmp == x

    @test Groebner._monom_overflow_check(ev)

    ev = Groebner.monom_construct_from_vector(PV{UInt8, Int, 5}, x)
    @test typeof(Groebner.monom_totaldeg(ev)) === UInt8
    @test Groebner.monom_totaldeg(ev) == UInt8(10)
    @test Groebner.monom_entrytype(ev) === Groebner.MonomHash

    for T in (UInt64, UInt32, UInt16, UInt8)
        for k in 1:10
            n = rand(1:20)
            x = rand(0:10, n)
            ev = Groebner.monom_construct_from_vector(PV{T, Int, n}, x)
            tmp = similar(x)
            @test x == Groebner.monom_to_vector!(tmp, ev)
            @test x == tmp
            @test Groebner.monom_entrytype(ev) === Groebner.MonomHash
        end
    end
end
