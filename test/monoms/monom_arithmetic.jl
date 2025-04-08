using Test, Random, Groebner

monom_construct_from_vector = Groebner.monom_construct_from_vector
monom_to_vector! = Groebner.monom_to_vector!
monom_copy = Groebner.monom_copy
monom_product! = Groebner.monom_product!
monom_division! = Groebner.monom_division!
monom_lcm! = Groebner.monom_lcm!
monom_is_gcd_const = Groebner.monom_is_gcd_const
monom_is_divisible = Groebner.monom_is_divisible
monom_is_divisible! = Groebner.monom_is_divisible!
monom_is_equal = Groebner.monom_is_equal

degree_types_to_test = [UInt64]
implementations_to_test = [
    Groebner.ExponentVector{T} where {T},
    Groebner.ExponentVector{UInt8},
    Groebner.PackedTuple1{T, UInt8} where {T},
    Groebner.PackedTuple2{T, UInt8} where {T},
    Groebner.PackedTuple3{T, UInt8} where {T},
    Groebner.PackedTuple4{T, UInt8} where {T}
]

@testset "monom arithmetic" begin
    for T in degree_types_to_test
        for EV in implementations_to_test
            MonomType = isconcretetype(EV) ? EV : EV{T}
            (T != UInt64) && (MonomType <: Groebner.AbstractPackedTuple) && continue

            a = monom_construct_from_vector(MonomType, [0, 0])
            b = monom_construct_from_vector(MonomType, [0, 0])
            c = monom_copy(a)
            @test monom_is_gcd_const(a, b)
            c = monom_product!(c, a, b)
            @test monom_is_equal(c, a)
            @test monom_is_equal(c, b)

            d = monom_construct_from_vector(MonomType, [1, 2])
            @test monom_is_equal(a, b)
            @test !monom_is_equal(a, d)
            @test monom_is_gcd_const(d, a)
            c = monom_product!(c, a, d)
            @test monom_is_equal(c, d)
            b = monom_product!(b, d, d)
            c = monom_product!(c, d, c)
            @test monom_is_equal(c, b)
            d = monom_product!(d, d, d)
            @test monom_is_equal(c, d)

            a = monom_construct_from_vector(MonomType, [2, 0])
            b = monom_construct_from_vector(MonomType, [0, 1])
            @test !monom_is_divisible(a, b)

            a = monom_construct_from_vector(MonomType, [0, 1, 3])
            b = monom_construct_from_vector(MonomType, [0, 0, 4])
            @test !monom_is_divisible(a, b)

            n = 6
            if n > Groebner.monom_max_vars(MonomType)
                continue
            end
            a = monom_construct_from_vector(MonomType, [1, 2, 0, 3, 0, 4])
            b = monom_construct_from_vector(MonomType, [0, 1, 0, 2, 0, 4])
            c = monom_construct_from_vector(MonomType, [1, 3, 0, 5, 0, 8])
            d = monom_construct_from_vector(MonomType, [1, 1, 0, 1, 0, 0])
            e = monom_construct_from_vector(MonomType, [1, 2, 0, 3, 0, 4])
            tmp = monom_copy(a)
            @test monom_is_divisible(a, b)
            flag, tmp = monom_is_divisible!(tmp, a, b)
            @test flag
            @test monom_is_equal(tmp, d)
            tmp = monom_product!(tmp, a, b)
            @test monom_is_equal(tmp, c)
            tmp = monom_division!(tmp, a, b)
            @test monom_is_equal(tmp, d)
            tmp = monom_lcm!(tmp, a, b)
            @test monom_is_equal(tmp, e)
            @test !monom_is_gcd_const(a, b)
            @test monom_is_equal(tmp, e)

            n = 9
            if n > Groebner.monom_max_vars(MonomType)
                continue
            end
            a = monom_construct_from_vector(MonomType, [1, 2, 0, 3, 0, 4, 0, 5, 6])
            b = monom_construct_from_vector(MonomType, [0, 1, 0, 2, 0, 4, 0, 0, 2])
            c = monom_construct_from_vector(MonomType, [1, 3, 0, 5, 0, 8, 0, 5, 8])
            d = monom_construct_from_vector(MonomType, [1, 1, 0, 1, 0, 0, 0, 5, 4])
            e = monom_construct_from_vector(MonomType, [1, 2, 0, 3, 0, 4, 0, 5, 6])
            tmp = monom_copy(a)
            @test monom_is_divisible(a, b)
            flag, tmp = monom_is_divisible!(tmp, a, b)
            @test flag
            @test monom_is_equal(tmp, d)
            tmp = monom_product!(tmp, a, b)
            @test monom_is_equal(tmp, c)
            tmp = monom_division!(tmp, a, b)
            @test monom_is_equal(tmp, d)
            tmp = monom_lcm!(tmp, a, b)
            @test monom_is_equal(tmp, e)
            @test !monom_is_gcd_const(a, b)
            @test monom_is_equal(tmp, e)

            n = 30
            if n > Groebner.monom_max_vars(MonomType)
                continue
            end
            a = monom_construct_from_vector(MonomType, 2 .* ones(Int, n))
            b = monom_construct_from_vector(MonomType, 1 .* ones(Int, n))
            tmp = monom_copy(a)
            @test monom_is_divisible(a, b)
            flag, tmp = monom_is_divisible!(tmp, a, b)
            @test flag
            @test monom_is_equal(tmp, b)
            @test !monom_is_gcd_const(a, b)
        end

        # test that different implementations agree
        for n in 1:5
            k = rand(1:50)

            implementations_to_test_spec = map(
                MonomType -> isconcretetype(MonomType) ? MonomType : MonomType{T},
                implementations_to_test
            )
            implementations_to_test_local = filter(
                MonomType -> Groebner.monom_max_vars(MonomType) >= k,
                implementations_to_test_spec
            )

            t = div(typemax(UInt8), 2k) - 1
            x, y = rand(1:t, k), rand(1:t, k)
            if sum(x) + sum(y) >= Groebner.monom_overflow_threshold(UInt8)
                continue
            end
            as = [monom_construct_from_vector(MT, x) for MT in implementations_to_test_local]
            bs = [monom_construct_from_vector(MT, y) for MT in implementations_to_test_local]

            results = []
            for (a, b) in zip(as, bs)
                c = monom_copy(a)
                tmp1 = similar(x)
                monom_to_vector!(tmp1, c)
                c = monom_product!(c, a, b)
                tmp2 = similar(x)
                monom_to_vector!(tmp2, c)
                flag = monom_is_divisible(a, b)
                _, a = monom_is_divisible!(a, c, b)
                tmp3 = similar(x)
                monom_to_vector!(tmp3, a)
                _, c = monom_is_divisible!(c, c, c)
                tmp4 = similar(x)
                monom_to_vector!(tmp4, a)
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
                if n > Groebner.monom_max_vars(MonomType)
                    continue
                end
                m = monom_construct_from_vector(MonomType, monom)
                ans = parse(Groebner.DivisionMask, mask, base=2)
                dm = Groebner.monom_create_divmask(
                    m,
                    Groebner.DivisionMask,
                    nvars,
                    divmap,
                    ndivbits,
                    false
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
                if n > Groebner.monom_max_vars(MonomType)
                    continue
                end
                t = div(typemax(UInt8), 8 * n)
                x, y = rand(1:max(t, 1), n), rand(1:max(t, 1), n)
                a = monom_construct_from_vector(MonomType, x)
                b = monom_construct_from_vector(MonomType, y)
                c = Groebner.monom_construct_const(MonomType, n)
                c = monom_product!(c, a, b)
                h = Groebner.monom_construct_hash_vector(Random.MersenneTwister(), MonomType, n)
                @test typeof(Groebner.monom_hash(a, h)) === Groebner.MonomHash
                @test Groebner.monom_hash(a, h) + Groebner.monom_hash(b, h) ==
                      Groebner.monom_hash(c, h)
            end
        end
    end
end
