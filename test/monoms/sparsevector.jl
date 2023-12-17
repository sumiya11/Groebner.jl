
@testset "sparse exponent vector" begin
    PV{T, N, I} = Groebner.SparseExponentVector{T, N, I} where {T, N, I}

    ev = Groebner.construct_const_monom(PV{UInt16, 60, Int}, 60)
    @test Groebner.totaldeg(ev) == UInt16(0)
    @test Groebner.monom_entrytype(ev) === Groebner.MonomHash

    x = [1, 2, 3, 0, 4]
    ev = Groebner.construct_monom(PV{UInt64, 5, Int}, x)
    @test typeof(Groebner.totaldeg(ev)) === UInt64
    @test Groebner.totaldeg(ev) == UInt64(10)
    @test Groebner.monom_entrytype(ev) === Groebner.MonomHash
    tmp = similar(x)
    @test Groebner.monom_to_dense_vector!(tmp, ev) == x
    @test tmp == x

    @test Groebner._monom_overflow_check(ev)

    ev = Groebner.construct_monom(PV{UInt8, 5, Int}, x)
    @test typeof(Groebner.totaldeg(ev)) === UInt8
    @test Groebner.totaldeg(ev) == UInt8(10)
    @test Groebner.monom_entrytype(ev) === Groebner.MonomHash

    for T in (UInt64, UInt32, UInt16, UInt8)
        for k in 1:10
            n = rand(1:20)
            x = rand(0:10, n)
            ev = Groebner.construct_monom(PV{T, n, Int}, x)
            tmp = similar(x)
            @test x == Groebner.monom_to_dense_vector!(tmp, ev)
            @test x == tmp
            @test Groebner.monom_entrytype(ev) === Groebner.MonomHash
        end
    end
end
