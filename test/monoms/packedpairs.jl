
@testset "packed exponent pair-1" begin
    EV{T, U} = Groebner.PackedPair1{T, U} where {T, U}

    x = [1, 2, 3, 0, 4]
    ev = Groebner.make_ev(EV{UInt64, UInt8}, x)
    @test ev.a1 == 0x0a00000400030201
    @test typeof(Groebner.totaldeg(ev)) === UInt64
    @test Groebner.totaldeg(ev) == UInt64(10)
    @test Groebner.powertype(ev) === Groebner.MonomHash
    tmp = similar(x)
    @test Groebner.make_dense!(tmp, ev) == [1, 2, 3, 0, 4]
    @test tmp == [1, 2, 3, 0, 4]

    @test Groebner.capacity(Groebner.PackedPair1{UInt64, UInt8}) == 7
    @test Groebner.capacity(Groebner.PackedPair1{UInt64, UInt16}) == 3

    y = [10, 20, 30, 40, 50]
    @test_throws Groebner.RecoverableException Groebner.make_ev(EV{UInt64, UInt8}, y)

    for B in (UInt8, UInt16,)
        for T in (UInt64, UInt32)
            for k in 1:10
                maxvars = div(sizeof(T), sizeof(B)) - 1
                n = rand(1:maxvars)
                @assert maxvars > 0
                x = rand(0:10, n)
                ev = Groebner.make_ev(EV{T, B}, x)
                tmp = similar(x)
                @test x == Groebner.make_dense!(tmp, ev)
                @test x == tmp
            end
        end
    end
end

@testset "packed exponent pair-2" begin
    EV{T, U} = Groebner.PackedPair2{T, U} where {T, U}

    x = [1, 2, 3, 0, 4]
    ev = Groebner.make_ev(EV{UInt64, UInt8}, x)
    @test ev.a1 == 0x0a00000400030201
    @test ev.a2 == 0x0000000000000000
    @test typeof(Groebner.totaldeg(ev)) === UInt64
    @test Groebner.totaldeg(ev) == UInt64(10)
    @test Groebner.powertype(ev) === Groebner.MonomHash
    tmp = similar(x)
    @test Groebner.make_dense!(tmp, ev) == [1, 2, 3, 0, 4]

    @test Groebner.capacity(Groebner.PackedPair2{UInt64, UInt8}) == 15
    @test Groebner.capacity(Groebner.PackedPair2{UInt64, UInt16}) == 7

    ev = Groebner.make_ev(EV{UInt64, UInt16}, x)
    @test (ev.a1, ev.a2) == (0x000a000000000004, 0x0000000300020001)
    @test typeof(Groebner.totaldeg(ev)) === UInt64
    @test Groebner.totaldeg(ev) == UInt16(10)
    
    y = [10, 20, 30, 40, 50]
    @test_throws Groebner.RecoverableException Groebner.make_ev(EV{UInt64, UInt8}, y)

    for B in (UInt8, UInt16,)
        for T in (UInt64, UInt32)
            for k in 1:10
                maxvars = 2*div(sizeof(T), sizeof(B)) - 1
                n = rand(1:maxvars)
                @assert maxvars > 0
                x = rand(0:10, n)
                tmp = similar(x)
                ev = Groebner.make_ev(EV{T, B}, x)
                @test x == Groebner.make_dense!(tmp, ev)
                @test x == tmp
            end
        end
    end
end

@testset "packed exponent pair-3" begin
    EV{T, U} = Groebner.PackedPair3{T, U} where {T, U}

    x = [1, 2, 3, 0, 4]
    ev = Groebner.make_ev(EV{UInt64, UInt8}, x)
    @test ev.a1 == 0x0a00000400030201
    @test ev.a2 == 0x0000000000000000
    @test typeof(Groebner.totaldeg(ev)) === UInt64
    @test Groebner.totaldeg(ev) == UInt64(10)
    @test Groebner.powertype(ev) === Groebner.MonomHash
    tmp = similar(x)
    @test Groebner.make_dense!(tmp, ev) == [1, 2, 3, 0, 4]

    @test Groebner.capacity(Groebner.PackedPair2{UInt64, UInt8}) == 15
    @test Groebner.capacity(Groebner.PackedPair2{UInt64, UInt16}) == 7

    ev = Groebner.make_ev(EV{UInt64, UInt16}, x)
    @test (ev.a1, ev.a2) == (0x000a000000000004, 0x0000000300020001)
    @test typeof(Groebner.totaldeg(ev)) === UInt64
    @test Groebner.totaldeg(ev) == UInt16(10)
    
    y = [10, 20, 30, 40, 50]
    @test_throws Groebner.RecoverableException Groebner.make_ev(EV{UInt64, UInt8}, y)

    y = [1, 2, 3, 4, 5]
    ev = Groebner.make_ev(EV{UInt64, UInt32}, y)
    @test ev.a1 == 0x0000000f00000005
    @test ev.a2 == 0x0000000400000003
    @test ev.a3 == 0x0000000200000001

    for B in (UInt8, UInt16,)
        for T in (UInt64, UInt32)
            for k in 1:10
                maxvars = 3*div(sizeof(T), sizeof(B)) - 1
                n = rand(1:maxvars)
                @assert maxvars > 0
                x = rand(0:div(typemax(B), 4*n), n)
                tmp = similar(x)
                ev = Groebner.make_ev(EV{T, B}, x)
                @test x == Groebner.make_dense!(tmp, ev)
                @test x == tmp
            end
        end
    end
end

