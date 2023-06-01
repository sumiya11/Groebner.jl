
make_ev = Groebner.make_ev
monom_product! = Groebner.monom_product!
monom_division! = Groebner.monom_division!
monom_lcm! = Groebner.monom_lcm!
is_gcd_const = Groebner.is_gcd_const
is_monom_divisible = Groebner.is_monom_divisible
is_monom_divisible! = Groebner.is_monom_divisible!
is_monom_elementwise_eq = Groebner.is_monom_elementwise_eq

B = UInt8
types_to_test = [UInt64, UInt32]
implementations_to_test = [
    Groebner.PowerVector{T} where {T},
    Groebner.PackedPair1{T, B} where {T},
    Groebner.PackedPair2{T, B} where {T},
    Groebner.PackedPair3{T, B} where {T},
    Groebner.PackedPair4{T, B} where {T}
    # Groebner.PackedVector{T, B} where {T},
]

@testset "term arithmetic" begin
    for T in types_to_test
        for EV in implementations_to_test
            a = make_ev(EV{T}, [0, 0])
            b = make_ev(EV{T}, [0, 0])
            c = copy(a)
            @test is_gcd_const(a, b)
            c = monom_product!(c, a, b)
            @test is_monom_elementwise_eq(c, a)
            @test is_monom_elementwise_eq(c, b)

            d = make_ev(EV{T}, [1, 2])
            @test is_monom_elementwise_eq(a, b)
            @test !is_monom_elementwise_eq(a, d)
            @test is_gcd_const(d, a)
            c = monom_product!(c, a, d)
            @test is_monom_elementwise_eq(c, d)
            b = monom_product!(b, d, d)
            c = monom_product!(c, d, c)
            @test is_monom_elementwise_eq(c, b)
            d = monom_product!(d, d, d)
            @test is_monom_elementwise_eq(c, d)

            n = 6
            if n > Groebner.capacity(EV{T})
                continue
            end
            a = make_ev(EV{T}, [1, 2, 0, 3, 0, 4])
            b = make_ev(EV{T}, [0, 1, 0, 2, 0, 4])
            c = make_ev(EV{T}, [1, 3, 0, 5, 0, 8])
            d = make_ev(EV{T}, [1, 1, 0, 1, 0, 0])
            e = make_ev(EV{T}, [1, 2, 0, 3, 0, 4])
            tmp = copy(a)
            @test is_monom_divisible(a, b)
            flag, tmp = is_monom_divisible!(tmp, a, b)
            @test flag
            @test is_monom_elementwise_eq(tmp, d)
            tmp = monom_product!(tmp, a, b)
            @test is_monom_elementwise_eq(tmp, c)
            tmp = monom_division!(tmp, a, b)
            @test is_monom_elementwise_eq(tmp, d)
            tmp = monom_lcm!(tmp, a, b)
            @test is_monom_elementwise_eq(tmp, e)
            @test !is_gcd_const(a, b)
            @test is_monom_elementwise_eq(tmp, e)

            n = 9
            if n > Groebner.capacity(EV{T})
                continue
            end
            a = make_ev(EV{T}, [1, 2, 0, 3, 0, 4, 0, 5, 6])
            b = make_ev(EV{T}, [0, 1, 0, 2, 0, 4, 0, 0, 2])
            c = make_ev(EV{T}, [1, 3, 0, 5, 0, 8, 0, 5, 8])
            d = make_ev(EV{T}, [1, 1, 0, 1, 0, 0, 0, 5, 4])
            e = make_ev(EV{T}, [1, 2, 0, 3, 0, 4, 0, 5, 6])
            tmp = copy(a)
            @test is_monom_divisible(a, b)
            flag, tmp = is_monom_divisible!(tmp, a, b)
            @test flag
            @test is_monom_elementwise_eq(tmp, d)
            tmp = monom_product!(tmp, a, b)
            @test is_monom_elementwise_eq(tmp, c)
            tmp = monom_division!(tmp, a, b)
            @test is_monom_elementwise_eq(tmp, d)
            tmp = monom_lcm!(tmp, a, b)
            @test is_monom_elementwise_eq(tmp, e)
            @test !is_gcd_const(a, b)
            @test is_monom_elementwise_eq(tmp, e)
        end
    end
end

@testset "term division mask" begin
    nvars = 4
    divmap = UInt32[
        0x00000001,
        0x00000002,
        0x00000003,
        0x00000004,
        0x00000005,
        0x00000006,
        0x00000007,
        0x00000008,
        0x00000001,
        0x00000002,
        0x00000003,
        0x00000004,
        0x00000005,
        0x00000006,
        0x00000007,
        0x00000008,
        0x00000001,
        0x00000002,
        0x00000003,
        0x00000004,
        0x00000005,
        0x00000006,
        0x00000007,
        0x00000008,
        0x00000001,
        0x00000002,
        0x00000003,
        0x00000004,
        0x00000005,
        0x00000006,
        0x00000007,
        0x00000008
    ]
    ndivbits = 8
    cases = [
        (monom=UInt8[0x00, 0x00, 0x00, 0x00], mask="00000000000000000000000000000000"),
        (monom=UInt8[0x01, 0x00, 0x1a, 0x00], mask="00000000111111110000000000000001"),
        (monom=UInt8[0x00, 0x01, 0x00, 0x00], mask="00000000000000000000000100000000"),
        (monom=UInt8[0x00, 0x00, 0x01, 0x00], mask="00000000000000010000000000000000"),
        (monom=UInt8[0x00, 0x00, 0x00, 0x01], mask="00000001000000000000000000000000"),
        (monom=UInt8[0x01, 0x01, 0x00, 0x00], mask="00000000000000000000000100000001"),
        (monom=UInt8[0x01, 0x00, 0x01, 0x00], mask="00000000000000010000000000000001"),
        (monom=UInt8[0x00, 0x03, 0x02, 0x00], mask="00000000000000110000011100000000"),
        (monom=UInt8[0x01, 0x00, 0x00, 0x01], mask="00000001000000000000000000000001"),
        (monom=UInt8[0x00, 0x01, 0x00, 0x01], mask="00000001000000000000000100000000"),
        (monom=UInt8[0x00, 0x00, 0x01, 0x01], mask="00000001000000010000000000000000"),
        (monom=UInt8[0x01, 0x01, 0x01, 0x00], mask="00000000000000010000000100000001"),
        (monom=UInt8[0x01, 0x01, 0x00, 0x01], mask="00000001000000000000000100000001"),
        (monom=UInt8[0x01, 0x0a, 0x01, 0x01], mask="00000001000000011111111100000001"),
        (monom=UInt8[0x00, 0x01, 0x01, 0x01], mask="00000001000000010000000100000000"),
        (monom=UInt8[0x01, 0x01, 0x01, 0x01], mask="00000001000000010000000100000001")
    ]

    for case in cases
        monom, mask = case.monom, case.mask
        for T in types_to_test
            for EV in implementations_to_test
                n = nvars
                if n > Groebner.capacity(EV{T})
                    continue
                end
                m = make_ev(EV{T}, monom)
                ans = parse(Groebner.DivisionMask, mask, base=2)
                dm = Groebner.monom_divmask(
                    m,
                    Groebner.DivisionMask,
                    nvars,
                    divmap,
                    ndivbits
                )
                @test ans == dm
            end
        end
    end
end

@testset "term hash linearity" begin
    for _ in 1:30
        for T in types_to_test
            for EV in implementations_to_test
                n = rand(1:50)
                if n > Groebner.capacity(EV{T})
                    continue
                end
                t = div(typemax(UInt8), 8 * n)
                x, y = rand(1:max(t, 1), n), rand(1:max(t, 1), n)
                a = make_ev(EV{T}, x)
                b = make_ev(EV{T}, y)
                c = Groebner.make_zero_ev(EV{T}, n)
                c = monom_product!(c, a, b)
                h = Groebner.make_hasher(EV{T}, n)
                @test typeof(hash(a, h)) === Groebner.MonomHash
                @test hash(a, h) + hash(b, h) == hash(c, h)
            end
        end
    end
end
