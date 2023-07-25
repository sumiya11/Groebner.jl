
construct_monom = Groebner.construct_monom
monom_to_dense_vector! = Groebner.monom_to_dense_vector!
copy_monom = Groebner.copy_monom
monom_product! = Groebner.monom_product!
monom_division! = Groebner.monom_division!
monom_lcm! = Groebner.monom_lcm!
is_gcd_const = Groebner.is_gcd_const
is_monom_divisible = Groebner.is_monom_divisible
is_monom_divisible! = Groebner.is_monom_divisible!
is_monom_elementwise_eq = Groebner.is_monom_elementwise_eq

degree_types_to_test = [UInt64, UInt32]
implementations_to_test = [
    Groebner.ExponentVector{T} where {T},
    Groebner.ExponentVector{UInt8},
    Groebner.PackedTuple1{T, UInt8} where {T},
    Groebner.PackedTuple2{T, UInt8} where {T},
    Groebner.PackedTuple3{T, UInt8} where {T},
    Groebner.SparseExponentVector{T} where {T}
]

@testset "monom arithmetic" begin
    for T in degree_types_to_test
        for EV in implementations_to_test
            MonomType = isconcretetype(EV) ? EV : EV{T}
            a = construct_monom(MonomType, [0, 0])
            b = construct_monom(MonomType, [0, 0])
            c = copy_monom(a)
            @test is_gcd_const(a, b)
            c = monom_product!(c, a, b)
            @test is_monom_elementwise_eq(c, a)
            @test is_monom_elementwise_eq(c, b)

            d = construct_monom(MonomType, [1, 2])
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

            a = construct_monom(MonomType, [2, 0])
            b = construct_monom(MonomType, [0, 1])
            @test !is_monom_divisible(a, b)

            a = construct_monom(MonomType, [0, 1, 3])
            b = construct_monom(MonomType, [0, 0, 4])
            @test !is_monom_divisible(a, b)

            n = 6
            if n > Groebner.max_vars_in_monom(MonomType)
                continue
            end
            a = construct_monom(MonomType, [1, 2, 0, 3, 0, 4])
            b = construct_monom(MonomType, [0, 1, 0, 2, 0, 4])
            c = construct_monom(MonomType, [1, 3, 0, 5, 0, 8])
            d = construct_monom(MonomType, [1, 1, 0, 1, 0, 0])
            e = construct_monom(MonomType, [1, 2, 0, 3, 0, 4])
            tmp = copy_monom(a)
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
            if n > Groebner.max_vars_in_monom(MonomType)
                continue
            end
            a = construct_monom(MonomType, [1, 2, 0, 3, 0, 4, 0, 5, 6])
            b = construct_monom(MonomType, [0, 1, 0, 2, 0, 4, 0, 0, 2])
            c = construct_monom(MonomType, [1, 3, 0, 5, 0, 8, 0, 5, 8])
            d = construct_monom(MonomType, [1, 1, 0, 1, 0, 0, 0, 5, 4])
            e = construct_monom(MonomType, [1, 2, 0, 3, 0, 4, 0, 5, 6])
            tmp = copy_monom(a)
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

        # test that different implementations agree
        for n in 1:5
            k = rand(1:50)

            implementations_to_test_spec = map(
                MonomType -> isconcretetype(MonomType) ? MonomType : MonomType{T},
                implementations_to_test
            )
            implementations_to_test_local = filter(
                MonomType -> Groebner.max_vars_in_monom(MonomType) >= k,
                implementations_to_test_spec
            )

            t = div(typemax(UInt8), 2k) - 1
            x, y = rand(1:t, k), rand(1:t, k)
            if sum(x) + sum(y) >= Groebner._monom_overflow_threshold(UInt8)
                continue
            end
            as = [construct_monom(MT, x) for MT in implementations_to_test_local]
            bs = [construct_monom(MT, y) for MT in implementations_to_test_local]

            results = []
            for (a, b) in zip(as, bs)
                c = copy_monom(a)
                tmp1 = similar(x)
                monom_to_dense_vector!(tmp1, c)
                c = monom_product!(c, a, b)
                tmp2 = similar(x)
                monom_to_dense_vector!(tmp2, c)
                flag = is_monom_divisible(a, b)
                _, a = is_monom_divisible!(a, c, b)
                tmp3 = similar(x)
                monom_to_dense_vector!(tmp3, a)
                _, c = is_monom_divisible!(c, c, c)
                tmp4 = similar(x)
                monom_to_dense_vector!(tmp4, a)
                push!(results, (tmp1, tmp2, flag, tmp3, tmp4))
            end
            @test length(unique(results)) == 1
        end
    end
end

@testset "monom division mask" begin
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
        for T in degree_types_to_test
            for EV in implementations_to_test
                MonomType = isconcretetype(EV) ? EV : EV{T}

                n = nvars
                if n > Groebner.max_vars_in_monom(MonomType)
                    continue
                end
                m = construct_monom(MonomType, monom)
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

@testset "monom hash linearity" begin
    for _ in 1:30
        for T in degree_types_to_test
            for EV in implementations_to_test
                MonomType = isconcretetype(EV) ? EV : EV{T}
                n = rand(1:50)
                if n > Groebner.max_vars_in_monom(MonomType)
                    continue
                end
                t = div(typemax(UInt8), 8 * n)
                x, y = rand(1:max(t, 1), n), rand(1:max(t, 1), n)
                a = construct_monom(MonomType, x)
                b = construct_monom(MonomType, y)
                c = Groebner.construct_const_monom(MonomType, n)
                c = monom_product!(c, a, b)
                h = Groebner.construct_hash_vector(MonomType, n)
                @test typeof(Groebner.monom_hash(a, h)) === Groebner.MonomHash
                @test Groebner.monom_hash(a, h) + Groebner.monom_hash(b, h) ==
                      Groebner.monom_hash(c, h)
            end
        end
    end
end
