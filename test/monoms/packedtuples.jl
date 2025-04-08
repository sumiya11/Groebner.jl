using Test, Groebner

@testset "packed exponent tuple-1" begin
    EV{T, U} = Groebner.PackedTuple1{T, U} where {T, U}

    x = [1, 2, 3, 0, 4]
    ev = Groebner.monom_construct_from_vector(EV{UInt64, UInt8}, x)
    @test ev.a1 == 0x0a00000400030201
    @test typeof(Groebner.monom_totaldeg(ev)) === UInt64
    @test Groebner.monom_totaldeg(ev) == UInt64(10)
    @test Groebner.monom_entrytype(ev) === UInt8
    tmp = similar(x)
    @test Groebner.monom_to_vector!(tmp, ev) == [1, 2, 3, 0, 4]
    @test tmp == [1, 2, 3, 0, 4]

    @test Groebner.monom_max_vars(Groebner.PackedTuple1{UInt64, UInt8}) == 7
    @test Groebner.monom_max_vars(Groebner.PackedTuple1{UInt64, UInt16}) == 3

    y = [10, 20, 30, 40, 50]
    @test_throws Groebner.MonomialDegreeOverflow Groebner.monom_construct_from_vector(
        EV{UInt64, UInt8},
        y
    )

    for B in (UInt8, UInt16)
        for T in (UInt64, UInt32)
            for k in 1:10
                maxvars = div(sizeof(T), sizeof(B)) - 1
                n = rand(1:maxvars)
                @assert maxvars > 0
                x = rand(0:10, n)
                ev = Groebner.monom_construct_from_vector(EV{T, B}, x)
                tmp = similar(x)
                @test x == Groebner.monom_to_vector!(tmp, ev)
                @test x == tmp
            end
        end
    end
end

@testset "packed exponent tuple-2" begin
    EV{T, U} = Groebner.PackedTuple2{T, U} where {T, U}

    x = [i for i in 1:15]
    ev = Groebner.monom_construct_from_vector(EV{UInt64, UInt8}, x)
    @test ev.a1 == 0x780f0e0d0c0b0a09
    @test ev.a2 == 0x0807060504030201
    @test Groebner.monom_totaldeg(ev) == sum(x)
    tmp = similar(x)
    @test Groebner.monom_to_vector!(tmp, ev) == x

    x = [1, 2, 3, 0, 4]
    ev = Groebner.monom_construct_from_vector(EV{UInt64, UInt8}, x)
    @test ev.a1 == 0x0a00000400030201
    @test ev.a2 == 0x0000000000000000
    @test typeof(Groebner.monom_totaldeg(ev)) === UInt64
    @test Groebner.monom_totaldeg(ev) == UInt64(10)
    @test Groebner.monom_entrytype(ev) === UInt8
    tmp = similar(x)
    @test Groebner.monom_to_vector!(tmp, ev) == [1, 2, 3, 0, 4]

    @test Groebner.monom_max_vars(Groebner.PackedTuple2{UInt64, UInt8}) == 15
    @test Groebner.monom_max_vars(Groebner.PackedTuple2{UInt64, UInt16}) == 7

    ev = Groebner.monom_construct_from_vector(EV{UInt64, UInt16}, x)
    @test (ev.a1, ev.a2) == (0x000a000000000004, 0x0000000300020001)
    @test typeof(Groebner.monom_totaldeg(ev)) === UInt64
    @test Groebner.monom_totaldeg(ev) == UInt16(10)

    y = [10, 20, 30, 40, 50]
    @test_throws Groebner.MonomialDegreeOverflow Groebner.monom_construct_from_vector(
        EV{UInt64, UInt8},
        y
    )

    for B in (UInt8, UInt16)
        for T in (UInt64, UInt32)
            for k in 1:10
                maxvars = 2 * div(sizeof(T), sizeof(B)) - 1
                n = rand(1:maxvars)
                @assert maxvars > 0
                x = rand(0:10, n)
                tmp = similar(x)
                ev = Groebner.monom_construct_from_vector(EV{T, B}, x)
                @test x == Groebner.monom_to_vector!(tmp, ev)
                @test x == tmp
            end
        end
    end
end

@testset "packed exponent tuple-3" begin
    EV{T, U} = Groebner.PackedTuple3{T, U} where {T, U}

    x = [div(i, 5) for i in 1:23]
    ev = Groebner.monom_construct_from_vector(EV{UInt64, UInt8}, x)
    @test ev.a1 == 0x2e04040404030303
    @test ev.a2 == 0x0303020202020201
    @test ev.a3 == 0x0101010100000000
    @test Groebner.monom_totaldeg(ev) == sum(x)
    tmp = similar(x)
    @test Groebner.monom_to_vector!(tmp, ev) == x

    x = [1, 2, 3, 0, 4]
    ev = Groebner.monom_construct_from_vector(EV{UInt64, UInt8}, x)
    @test ev.a1 == 0x0a00000400030201
    @test ev.a2 == 0x0000000000000000
    @test ev.a3 == 0x0000000000000000
    @test typeof(Groebner.monom_totaldeg(ev)) === UInt64
    @test Groebner.monom_totaldeg(ev) == UInt64(10)
    @test Groebner.monom_entrytype(ev) === UInt8
    tmp = similar(x)
    @test Groebner.monom_to_vector!(tmp, ev) == [1, 2, 3, 0, 4]

    @test Groebner.monom_max_vars(Groebner.PackedTuple3{UInt64, UInt8}) == 23
    @test Groebner.monom_max_vars(Groebner.PackedTuple3{UInt64, UInt16}) == 11

    ev = Groebner.monom_construct_from_vector(EV{UInt64, UInt16}, x)
    @test (ev.a1, ev.a2) == (0x000a000000000004, 0x0000000300020001)
    @test typeof(Groebner.monom_totaldeg(ev)) === UInt64
    @test Groebner.monom_totaldeg(ev) == UInt16(10)

    y = [10, 20, 30, 40, 50]
    @test_throws Groebner.MonomialDegreeOverflow Groebner.monom_construct_from_vector(
        EV{UInt64, UInt8},
        y
    )

    y = [1, 2, 3, 4, 5]
    ev = Groebner.monom_construct_from_vector(EV{UInt64, UInt32}, y)
    @test ev.a1 == 0x0000000f00000005
    @test ev.a2 == 0x0000000400000003
    @test ev.a3 == 0x0000000200000001

    for B in (UInt8, UInt16)
        for T in (UInt64, UInt32)
            for k in 1:10
                maxvars = 3 * div(sizeof(T), sizeof(B)) - 1
                n = rand(1:maxvars)
                @assert maxvars > 0
                x = rand(0:div(typemax(B), 4 * n), n)
                tmp = similar(x)
                ev = Groebner.monom_construct_from_vector(EV{T, B}, x)
                @test x == Groebner.monom_to_vector!(tmp, ev)
                @test x == tmp
            end
        end
    end
end

@testset "packed exponent tuple-4" begin
    EV{T, U} = Groebner.PackedTuple4{T, U} where {T, U}

    x = [1, 2, 3, 0, 4]
    ev = Groebner.monom_construct_from_vector(EV{UInt64, UInt8}, x)
    @test ev.a1 == 0x0a00000400030201
    @test ev.a2 == 0x0000000000000000
    @test ev.a3 == 0x0000000000000000
    @test ev.a4 == 0x0000000000000000
    @test typeof(Groebner.monom_totaldeg(ev)) === UInt64
    @test Groebner.monom_totaldeg(ev) == UInt64(10)
    @test Groebner.monom_entrytype(ev) === UInt8
    tmp = similar(x)
    @test Groebner.monom_to_vector!(tmp, ev) == [1, 2, 3, 0, 4]

    @test Groebner.monom_max_vars(Groebner.PackedTuple4{UInt64, UInt8}) == 31
    @test typeof(Groebner.monom_totaldeg(ev)) === UInt64
    @test Groebner.monom_totaldeg(ev) == UInt16(10)

    y = [10, 20, 30, 40, 50]
    @test_throws Groebner.MonomialDegreeOverflow Groebner.monom_construct_from_vector(
        EV{UInt64, UInt8},
        y
    )

    y = [1, 2, 3, 4, 5]
    ev = Groebner.monom_construct_from_vector(EV{UInt64, UInt32}, y)
    @test ev.a1 == 0x0000000f00000005
    @test ev.a2 == 0x0000000400000003
    @test ev.a3 == 0x0000000200000001

    for B in (UInt8, UInt16)
        for T in (UInt64, UInt32)
            for k in 1:10
                maxvars = 4 * div(sizeof(T), sizeof(B)) - 1
                n = rand(1:maxvars)
                @assert maxvars > 0
                x = rand(0:div(typemax(B), 4 * n), n)
                tmp = similar(x)
                ev = Groebner.monom_construct_from_vector(EV{T, B}, x)
                @test x == Groebner.monom_to_vector!(tmp, ev)
                @test x == tmp
            end
        end
    end
end
