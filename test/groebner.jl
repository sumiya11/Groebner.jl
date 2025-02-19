using AbstractAlgebra
using Base.Threads
using Combinatorics, Primes, Random
using Test

@testset "groebner basic" begin
    R, (x, y) = polynomial_ring(GF(2^31 - 1), ["x", "y"], internal_ordering=:degrevlex)

    @test_throws DomainError Groebner.groebner([])
    @test_throws DomainError Groebner.groebner([1])
    @test_throws DomainError Groebner.groebner(Vector{elem_type(R)}())

    @test Groebner.groebner([R(0)]) == [R(0)]
    @test Groebner.groebner([x, x, x]) == [x]
    @test Groebner.groebner([x, x, y, y, y, x, x, y]) == [y, x]
    @test Groebner.groebner([R(1)]) == [R(1)]
    @test Groebner.groebner([x^2 + y, y * x - 1, R(1), y^4]) == [R(1)]
    @test Groebner.groebner([2x, 3y, 4x + 5y]) == [y, x]

    R, (x1, x2) = polynomial_ring(GF(2^31 - 1), ["x1", "x2"], internal_ordering=:degrevlex)

    fs = [x1^2 * x2^2 + 2 * x1^2 * x2, x1^2 * x2^2 + 3 * x1^2 * x2 + 5 * x1 * x2^2, x1^2 * x2^2]
    @test Groebner.groebner(fs) == [x1 * x2^2, x1^2 * x2]

    fs = [x1 * x2^2 + x1 * x2, x1^2 * x2^2 + 2 * x1 * x2, 2 * x1^2 * x2^2 + x1 * x2]
    @test Groebner.groebner(fs) == [x1 * x2]

    @test_throws AssertionError Groebner.groebner([x1, x2], this_keyword_does_not_exist=":^(")
end

get_data(sys, T) =
    (map(f -> collect(exponent_vectors(f)), sys), map(f -> collect(T.(coefficients(f))), sys))

function test_low_level_interface(ring, sys; passthrough...)
    T(x) = base_ring(parent(sys[1])) == QQ ? Rational{BigInt}(x) : UInt64(lift(x))
    monoms, coeffs = get_data(sys, T)
    gb_monoms1, gb_coeffs1 = Groebner.groebner(ring, monoms, coeffs; passthrough...)
    gb = Groebner.groebner(sys; passthrough...)
    gb_monoms2, gb_coeffs2 = get_data(gb, T)
    @test (gb_monoms1, gb_coeffs1) == (gb_monoms2, gb_coeffs2)
    @test Groebner.isgroebner(ring, gb_monoms1, gb_coeffs1)
    @test all(iszero, Groebner.normalform(ring, gb_monoms1, gb_coeffs1, ring, monoms, coeffs))
end

@testset "groebner low level" begin
    R, (x, y) = polynomial_ring(GF(2^31 - 1), ["x", "y"], internal_ordering=:degrevlex)
    ring_ff = Groebner.PolyRing(2, Groebner.DegRevLex(), 2^31 - 1)
    ring_ff2 = Groebner.PolyRing(2, Groebner.DegRevLex(), 2^32 - 5)
    ring_qq = Groebner.PolyRing(2, Groebner.DegRevLex(), 0)

    @test_throws DomainError Groebner.groebner(
        Groebner.PolyRing(0, Groebner.DegRevLex(), 0),
        [[[0, 0]]],
        [[1]]
    )
    @test_throws DomainError Groebner.groebner(
        Groebner.PolyRing(2, Groebner.DegRevLex(), -2),
        [[[0, 0]]],
        [[1]]
    )
    @test_throws DomainError Groebner.groebner(
        Groebner.PolyRing(2, Groebner.InputOrdering(), 2^31 - 1),
        [[[0, 0]]],
        [[1]]
    )
    @test_throws DomainError Groebner.groebner(
        Groebner.PolyRing(2, Groebner.DegRevLex(), 33),
        [[[0, 0]]],
        [[1]]
    )
    @test_throws DomainError Groebner.groebner(
        Groebner.PolyRing(2, Groebner.DegRevLex(3, 2, 1), 0),
        [[[0, 0]]],
        [[1]]
    )
    @test_throws DomainError Groebner.groebner(
        Groebner.PolyRing(2, Groebner.DegRevLex(1, 2, 2), 0),
        [[[0, 0]]],
        [[1]]
    )

    @test_throws DomainError Groebner.groebner(ring_ff, [], [])
    @test_throws DomainError Groebner.groebner(ring_ff, [[]], [[]])
    @test_throws DomainError Groebner.groebner(ring_ff, [[[1, 2], []]], [[1]])
    @test_throws DomainError Groebner.groebner(ring_ff, [[[1, 2, 3]]], [[1]])
    @test_throws DomainError Groebner.groebner(ring_ff, [[[1, 2]]], [[1, 2]])
    @test_throws DomainError Groebner.groebner(ring_ff, [[[1, 2]]], [[1, 2]])

    @test ([Vector{Vector{Int}}()], [Int[]]) ==
          Groebner.groebner(ring_ff, [Vector{Vector{Int}}(), Vector{Vector{Int}}()], [Int[], Int[]])
    @test ([[[0, 0]]], [[1]]) == Groebner.groebner(ring_ff, [[[0, 0]]], [[1]])
    @test ([[[1, 1]]], [[1]]) == Groebner.groebner(ring_ff, [[[1, 1]]], [[2]])
    @test ([[[1, 1]]], [[1]]) ==
          Groebner.groebner(ring_ff, [[[1, 1]]], [[2]]; ordering=Groebner.DegRevLex())
    @test ([[[1, 1]]], [[1]]) ==
          Groebner.groebner(ring_ff, [[[1, 1]]], [[2]]; ordering=Groebner.Lex(1, 2))
    @test ([[[0, 1], [0, 0]]], [[1, 5]]) ==
          Groebner.groebner(ring_ff, [[[0, 1], [0, 0], [0, 0]]], [[1, 2, 3]])
    @test ([[[0, 1], [0, 0]]], [[1, 5]]) ==
          Groebner.groebner(ring_qq, [[[0, 1], [0, 0], [0, 0]]], [[1, 2, 3]])
    @test ([[[0, 1], [0, 0]]], [[1, 2^31 - 2]]) ==
          Groebner.groebner(ring_ff, [[[0, 1], [0, 0], [0, 0]]], [[1, 2, -3]])
    @test ([[[0, 1]]], [[1]]) ==
          Groebner.groebner(ring_ff, [[[0, 1], [0, 0], [0, 0]]], [[1, 3, -3]])
    @test ([Vector{Vector{Int}}()], [Int[]]) == Groebner.groebner(ring_ff, [[[0, 1]]], [[0]])
    @test ([Vector{Vector{Int}}()], [Int[]]) == Groebner.groebner(ring_qq, [[[0, 1]]], [[0]])
    @test ([[[1, 1]]], [[1]]) == Groebner.groebner(ring_ff, [[[1, 1]]], [[2^31]])
    @test ([[[1, 1], [1, 0]]], [[1, 5]]) == Groebner.groebner(
        ring_ff2,
        [[[1, 1], [1, 0], [1, 0], [0, 0], [0, 0], [0, 0]]],
        [[1, 2^31, 2^31, 2^31, 2^31, -5]]
    )

    # Sorting on the fly
    ring_qq = Groebner.PolyRing(2, Groebner.DegRevLex(), 0)
    @test ([[[1, 1], [0, 0]]], [[1, 3 // 4]]) ==
          Groebner.groebner(ring_qq, [[[0, 0], [1, 1]]], [[3, 4]])

    ring_gf = Groebner.PolyRing(5, Groebner.Lex(), 2^30 + 3)
    sys = Groebner.Examples.cyclicn(5, k=GF(2^30 + 3), internal_ordering=:lex)
    perm = map(i -> shuffle(collect(1:length(sys[i]))), 1:length(sys))
    monoms = map(i -> collect(exponent_vectors(sys[i]))[perm[i]], 1:length(sys))
    coeffs = map(i -> collect(UInt64.(lift.(coefficients(sys[i]))))[perm[i]], 1:length(sys))
    @test get_data(Groebner.groebner(sys), UInt64 ∘ lift) ==
          Groebner.groebner(ring_gf, monoms, coeffs)

    ring = Groebner.PolyRing(2, Groebner.Lex(1, 2), 0)
    @test ([[[1, 0], [0, 1]]], [[1, 6 // 5]]) ==
          Groebner.groebner(ring, [[[1, 0], [0, 1]]], [[-5, -6]])
    ring = Groebner.PolyRing(2, Groebner.Lex(2, 1), 0)
    @test ([[[0, 1], [1, 0]]], [[1, 5 // 6]]) ==
          Groebner.groebner(ring, [[[1, 0], [0, 1]]], [[-5, -6]])
    ring = Groebner.PolyRing(2, Groebner.DegRevLex(1, 2), 7)
    @test ([[[2, 2], [0, 3]]], [[1, 4]]) ==
          Groebner.groebner(ring, [[[0, 3], [2, 2], [2, 1]]], [[1, 2, 0]])

    sys = Groebner.Examples.cyclicn(5, k=GF(2^40 + 15))
    ring = Groebner.PolyRing(5, Groebner.DegRevLex(), 2^40 + 15)
    test_low_level_interface(ring, sys)

    sys = Groebner.Examples.cyclicn(5, k=GF(2^30 + 3), internal_ordering=:lex)
    ring = Groebner.PolyRing(5, Groebner.Lex(), 2^30 + 3)
    test_low_level_interface(ring, sys)

    sys = Groebner.Examples.cyclicn(5, k=QQ)
    ring = Groebner.PolyRing(5, Groebner.DegRevLex(), 0)
    test_low_level_interface(ring, sys)

    sys = Groebner.Examples.cyclicn(5, k=QQ)
    sys_lex = Groebner.Examples.cyclicn(5, k=QQ; internal_ordering=:lex)
    ring = Groebner.PolyRing(5, Groebner.Lex(), 0)
    monoms, coeffs = get_data(sys, Rational)
    gb = get_data(Groebner.groebner(sys_lex), Rational)
    @test gb == Groebner.groebner(ring, monoms, coeffs)

    sys = Groebner.Examples.cyclicn(7, k=GF(11); internal_ordering=:lex)
    sys_drl = Groebner.Examples.cyclicn(7, k=GF(11); internal_ordering=:degrevlex)
    ring = Groebner.PolyRing(7, Groebner.DegRevLex(), 11)
    monoms, coeffs = get_data(sys, UInt64 ∘ lift)
    gb = get_data(Groebner.groebner(sys_drl), UInt64 ∘ lift)
    @test gb == Groebner.groebner(ring, monoms, coeffs)
end

@testset "groebner generic" begin
    # GB over
    # - Zp
    # - Zp for large p
    # - QQ
    n = 5
    syss = [
        Groebner.Examples.cyclicn(n, k=GF(2^40 + 15)),
        Groebner.Examples.cyclicn(n, k=GF(2^40 + 15), internal_ordering=:lex),
        Groebner.Examples.cyclicn(n, k=GF(2)),
        Groebner.Examples.cyclicn(n, k=GF(nextprime(big(2)^1000))),
        Groebner.Examples.cyclicn(n, k=GF(nextprime(big(2)^1000)), internal_ordering=:lex),
        Groebner.Examples.cyclicn(n, k=QQ)
    ]
    for sys in syss
        ord = internal_ordering(parent(sys[1])) == :lex ? Groebner.Lex() : Groebner.DegRevLex()
        ring = Groebner.PolyRing(n, ord, 0, :generic)
        exps, cfs = get_data(sys, Groebner.CoeffGeneric)
        gbexps, gbcfs = groebner(ring, exps, cfs)
        gb = groebner(sys)
        @test isgroebner(gb)
        @test gbexps == get_data(gb, identity)[1]
        @test isgroebner(ring, gbexps, gbcfs)
        @test all(isempty, normalform(ring, gbexps, gbcfs, ring, exps, cfs)[1])
        @test map(c -> Groebner.generic_coeff_data.(c), gbcfs) == get_data(gb, identity)[2]
    end

    # GB over QQ(a,b)
    R_param, (a, b) = QQ["a", "b"]
    R, (x, y, z) = fraction_field(R_param)["x", "y", "z"]
    F = [x^2 + x + (a + 1), x * y + b * y * z + 1 // (a * b), x * z + z + b]
    @test groebner(F) == [
        z^2 + b // (a + 1) * z + b^2 // (a + 1),
        y + (-a - 1) // (a^2 * b^2 + a * b^4 + a * b^2) * z - 1 // (a^2 * b + a * b^3 + a * b),
        x + (-a - 1) // b * z
    ]

    @test groebner([R(a - b), R(a)]) == [one(R)]
end

@testset "groebner reduced=false" begin
    R, (x, y) = polynomial_ring(GF(2^31 - 1), ["x", "y"], internal_ordering=:lex)

    fs = [x + y^2, x * y - y^2]
    gb = Groebner.groebner(fs, reduced=false)
    @test Groebner.isgroebner(gb)

    fs = [x + y, x^2 + y]
    gb = Groebner.groebner(fs, reduced=false)
    @test Groebner.isgroebner(gb)

    fs = [y, x * y + x]
    gb = Groebner.groebner(fs, reduced=false)
    @test gb == [y, x]

    fs = [x^2 + 5, 2y^2 + 3]
    gb = Groebner.groebner(fs, reduced=false)
    @test gb == [y^2 + 1073741825, x^2 + 5]

    fs = [y^2 + x, x^2 * y + y, y + 1]
    gb = Groebner.groebner(fs, reduced=false)
    @test gb == [R(1)]

    fs = [x^2 * y^2, 2 * x * y^2 + 3 * x * y]
    gb = Groebner.groebner(fs, reduced=false)
    @test gb == [x * y^2 + 1073741825 * x * y, x^2 * y]
end

@testset "groebner ground fields" begin
    # Large fields
    fields = [
        GF(2^20 + 7),
        GF(2^29 - 3),
        GF(2^29 + 11),
        GF(2^30 - 35),
        GF(2^30 + 3),
        GF(2^31 - 1),
        GF(2^31 + 11),
        GF(2^32 + 15)
    ]
    fields = vcat(fields, map(GF, Primes.nextprimes(2^20, 10)))
    for field in fields
        R, (x, y) = field["x", "y"]

        @test Groebner.groebner([R(3), R(4)]) == [R(1)]
        @test Groebner.groebner([x - 3, y + 3]) == [y + 3, x - 3]
        noon = Groebner.Examples.noonn(5, k=field, internal_ordering=:degrevlex)
        gb = Groebner.groebner(noon)
        @test Groebner.isgroebner(gb) && all(iszero, Groebner.normalform(gb, noon))

        trace, gb00 = Groebner.groebner_learn(gb)
        flag, gb01 = Groebner.groebner_apply!(trace, gb)
        @test gb == gb00 == gb01 && flag
    end

    # Larger fields
    fields =
        [GF(BigInt(2)^40 + 15), GF(BigInt(2)^63 - 25), GF(BigInt(2)^63 + 29), GF(BigInt(2)^64 - 59)]
    for field in fields
        R, (x, y) = field["x", "y"]

        @test Groebner.groebner([R(3), R(4)]) == [R(1)]
        @test Groebner.groebner([x - 3, y + 3]) == [y + 3, x - 3]
        noon = Groebner.Examples.noonn(5, k=field, internal_ordering=:degrevlex)
        gb = Groebner.groebner(noon)
        @test Groebner.isgroebner(gb) && all(iszero, Groebner.normalform(gb, noon))

        trace, gb00 = Groebner.groebner_learn(gb)
        flag, gb01 = Groebner.groebner_apply!(trace, gb)
        @test gb == gb00 == gb01 && flag
    end

    # Smaller fields
    fields = map(GF, [113, 127, 131, 251, 257, 32749, 32771, 65521, 65537])
    for field in fields
        R, (x, y) = field["x", "y"]

        @test Groebner.groebner([R(3), R(4)]) == [R(1)]
        @test Groebner.groebner([x - 3, y + 3]) == [y + 3, x - 3]
        noon = Groebner.Examples.noonn(5, k=field, internal_ordering=:degrevlex)
        gb = Groebner.groebner(noon)
        @test Groebner.isgroebner(gb) && all(iszero, Groebner.normalform(gb, noon))

        trace, gb00 = Groebner.groebner_learn(gb)
        flag, gb01 = Groebner.groebner_apply!(trace, gb)
        @test gb == gb00 == gb01 && flag
    end

    # Tiny fields
    fields = map(GF, [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41])
    for field in fields
        R, (x, y) = field["x", "y"]

        @test Groebner.groebner([R(3), R(4)]) == [R(1)]
        @test Groebner.groebner([x - 3, y + 3]) == [y + 3, x - 3]
        noon = Groebner.Examples.noonn(5, k=field, internal_ordering=:degrevlex)
        gb = Groebner.groebner(noon)
        @test Groebner.isgroebner(gb) && all(iszero, Groebner.normalform(gb, noon))

        trace, gb00 = Groebner.groebner_learn(gb)
        flag, gb01 = Groebner.groebner_apply!(trace, gb)
        @test gb == gb00 == gb01 && flag
    end
end

@testset "groebner modular" begin
    for modular in [:classic_modular, :learn_and_apply]
        R, (x,) = polynomial_ring(QQ, ["x"], internal_ordering=:degrevlex)
        @test Groebner.groebner([x], modular=modular) == [x]
        @test Groebner.groebner([4x], modular=modular) == [x]
        @test Groebner.groebner([x + 1], modular=modular) == [x + 1]
        @test Groebner.groebner([4x + 1], modular=modular) == [x + 1 // 4]
        @test Groebner.groebner([x + 5], modular=modular) == [x + 5]
        @test Groebner.groebner([x + 2^10], modular=modular) == [x + 2^10]
        @test Groebner.groebner([x + 2^30], modular=modular) == [x + 2^30]
        @test Groebner.groebner([x + 2^30 + 3], modular=modular) == [x + 2^30 + 3]
        @test Groebner.groebner([x + 2^31 - 1], modular=modular) == [x + 2^31 - 1]
        @test Groebner.groebner([x + 2^31 - 1, x^2], modular=modular) == [1]
        @test Groebner.groebner([(3232323 // 7777)x + 7777 // 3232323], modular=modular) ==
              [x + 60481729 // 10447911976329]
        @test Groebner.groebner([((2^31 - 1) // 1)x + 1], modular=modular) == [x + 1 // 2147483647]
        @test Groebner.groebner([(1 // (2^31 - 1))x + 1], modular=modular) == [x + 2147483647]
        @test Groebner.groebner(
            [1 // (2^30 + 3) * x^2 + (2^30 + 3)x + 1 // (1073741831)],
            modular=modular
        ) == [x^2 + 1152921511049297929x + (2^30 + 3) // 1073741831]

        R, (x, y) = polynomial_ring(QQ, ["x", "y"], internal_ordering=:degrevlex)

        fs = [x, y]
        G = Groebner.groebner(fs, modular=modular)
        @test G == [y, x]

        fs = [5y, 3 // 4 * x]
        G = Groebner.groebner(fs, modular=modular)
        @test G == [y, x]

        fs = [x^2 - (1 // 6) * y, x * y]
        G = Groebner.groebner(fs, modular=modular)
        @test G == [y^2, x * y, x^2 - (1 // 6)y]

        fs = [
            QQ(11, 3) * x^2 * y - QQ(2, 4) * x * y - QQ(1, 7) * y,
            QQ(1, 1) * x * y + QQ(7, 13) * x
        ]
        G = Groebner.groebner(fs, modular=modular)
        @test G == [y^2 + 7 // 13 * y, x * y + 7 // 13 * x, x^2 - 3 // 22 * x + 39 // 539 * y]

        root = Groebner.Examples.rootn(3, k=QQ, internal_ordering=:degrevlex)
        gb = Groebner.groebner(root, modular=modular)
        @test Groebner.isgroebner(gb)

        root = Groebner.Examples.rootn(4, k=QQ, internal_ordering=:degrevlex)
        gb = Groebner.groebner(root, modular=modular)
        @test Groebner.isgroebner(gb)

        root = Groebner.Examples.rootn(8, k=QQ, internal_ordering=:degrevlex)
        gb = Groebner.groebner(root, modular=modular)
        @test Groebner.isgroebner(gb)

        noon = Groebner.Examples.noonn(3, k=QQ, internal_ordering=:degrevlex)
        gb = Groebner.groebner(noon, modular=modular)
        @test Groebner.isgroebner(gb)

        noon = Groebner.Examples.noonn(4, k=QQ, internal_ordering=:degrevlex)
        gb = Groebner.groebner(noon, modular=modular)
        @test Groebner.isgroebner(gb)

        chan = Groebner.Examples.chandran(6)
        gb_truth = Groebner.groebner(chan)
        for _composite in (1, 2, 4, 8, 16)
            gb = groebner(chan, _composite=_composite)
            @test gb == gb_truth
        end

        # Test a number of cases directly
        R, (x, y, z) = QQ["x", "y", "z"]
        xs = gens(R)
        system = [z + 1, -y^5 + y^2 + y + 1, x^10 + x^9 + x^7 - x^5 + x^2 + x]
        for i in 1:100
            # Substitute a point of hight of up to 1000 bits. When raised to power
            # 10, this becomes at most 10k bits
            point = map(x -> x * BigInt(2)^rand(1:(10i)), xs)
            system_i = map(poly -> evaluate(poly, point), system)
            gb_i = Groebner.groebner(system_i, modular=modular)
            @test gb_i == map(poly -> divexact(poly, leading_coefficient(poly)), system_i)
        end
    end

    R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], internal_ordering=:degrevlex)

    fs = [(12345678 // 12347)x, (222222221111123 // 2131232232097)y + z]
    G = Groebner.groebner(fs)
    @test G == [y + 2131232232097 // 222222221111123 * z, x]

    fs = [(224 // 225) * x + y, x^2 + 1]
    ans = [x + 225 // 224 * y, y^2 + 50176 // 50625]
    @test Groebner.groebner(fs) == ans

    fs = [(12345 // 12345678)x + y, x^2 + 1]
    ans = [x + 4115226 // 4115 * y, y^2 + 16933225 // 16935085031076]
    @test Groebner.groebner(fs) == ans

    fs = [(12345 // 12345678)x + y, x^2 + z^2, y^2 + w^2, x + 2 * y + 3 * z + 4 * w]
    G = Groebner.groebner(fs)

    @test G == [
        y - 12345 // 4106996 * z - 4115 // 1026749 * w,
        x + 6172839 // 2053498 * z + 4115226 // 1026749 * w,
        z * w + 1692834523553 // 4063974 * w^2,
        z^2 - 16935085031076 // 16933225 * w^2,
        w^3
    ]

    G = [(2^30 + 3)x, (1 // (2^30 + 3))y]
    G = Groebner.groebner(G)
    @test G == [y, x]

    G = [(2^31 - 1) * x + y, (1 // (2^31 - 1)) * y + 1]
    G = Groebner.groebner(G)
    @test G == [y + 2147483647, x - 1]

    R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], internal_ordering=:lex)
    fs = [(12345678 // 12347)x, (222222221111123 // 2131232232097)y + z]
    G = Groebner.groebner(fs)
    @test G == [y + 2131232232097 // 222222221111123 * z, x]

    P = prod(BigInt, prevprimes(2^31, 100))
    @test groebner([x^2 + P * x + 1]) == [x^2 + P * x + 1]
end

@testset "groebner output sorted" begin
    for K in [GF(11), GF(307), GF(12007), GF(2^31 - 1), QQ]
        R, (x, y, z) = polynomial_ring(K, ["x", "y", "z"], internal_ordering=:lex)

        @test Groebner.groebner([x, y]) == Groebner.groebner([y, x]) == [y, x]
        @test Groebner.groebner([x^2 + y, x * y]) == [y^2, x * y, x^2 + y]
        @test Groebner.groebner([3x + 2, 5y]) == [y, x + K(2) // K(3)]
    end
end

@testset "monomial overflow" begin
    R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"], internal_ordering=:degrevlex)

    @test groebner([x^(2^31)]) == [x^2^31]
    @test_throws Groebner.MonomialDegreeOverflow groebner([x^(2^33)])

    for monoms in [:auto, :dense, :packed]
        gb_1 = [x * y^100 + y, x^100 * y + y^100, y^199 + 2147483646 * x^99 * y]
        gb_2 = [x * y^200 + y, x^200 * y + y^200, y^399 + 2147483646 * x^199 * y]
        gb_3 = [x * y^1000 + y, x^1000 * y + y^1000, y^1999 + 2147483646 * x^999 * y]
        @test Groebner.groebner([x^100 * y + y^100, x * y^100 + y], monoms=monoms) == gb_1
        @test Groebner.groebner([x^200 * y + y^200, x * y^200 + y], monoms=monoms) == gb_2
        @test Groebner.groebner([x^1000 * y + y^1000, x * y^1000 + y], monoms=monoms) == gb_3

        @test Groebner.isgroebner(gb_1)
        @test Groebner.isgroebner(gb_2)
        @test Groebner.isgroebner(gb_3)

        @test Groebner.normalform(gb_1, [x, y, R(1), R(0), x^1000]) == [x, y, R(1), R(0), x^1000]
        @test Groebner.normalform(gb_2, [x, y, R(1), R(0), x^1000]) == [x, y, R(1), R(0), x^1000]
        @test Groebner.normalform(gb_3, [x, y, R(1), R(0), x^10]) == [x, y, R(1), R(0), x^10]
        @test Groebner.normalform(gb_3, [x, y, R(1), R(0), x^1000]) == [x, y, R(1), R(0), x^1000]
    end
end

@testset "groebner reduced=true" begin
    root = Groebner.Examples.rootn(3, k=GF(2^31 - 1), internal_ordering=:degrevlex)
    x1, x2, x3 = gens(parent(first(root)))
    gb = Groebner.groebner(root, reduced=true)
    @test gb == [x1 + x2 + x3, x2^2 + x2 * x3 + x3^2, x3^3 + 2147483646]

    root = Groebner.Examples.rootn(6, k=GF(2^31 - 1), internal_ordering=:degrevlex)
    x1, x2, x3, x4, x5, x6 = gens(parent(first(root)))
    gb = Groebner.groebner(root, reduced=true)
    #! format: off
    @test gb == [
        x1 + x2 + x3 + x4 + x5 + x6,
        x2^2 + x2*x3 + x3^2 + x2*x4 + x3*x4 + x4^2 + x2*x5 + x3*x5 + x4*x5 + x5^2 + x2*x6 + x3*x6 + x4*x6 + x5*x6 + x6^2,
        x3^3 + x3^2*x4 + x3*x4^2 + x4^3 + x3^2*x5 + x3*x4*x5 + x4^2*x5 + x3*x5^2 + x4*x5^2 + x5^3 + x3^2*x6 + x3*x4*x6 + x4^2*x6 + x3*x5*x6 + x4*x5*x6 + x5^2*x6 + x3*x6^2 + x4*x6^2 + x5*x6^2 + x6^3,
        x4^4 + x4^3*x5 + x4^2*x5^2 + x4*x5^3 + x5^4 + x4^3*x6 + x4^2*x5*x6 + x4*x5^2*x6 + x5^3*x6 + x4^2*x6^2 + x4*x5*x6^2 + x5^2*x6^2 + x4*x6^3 + x5*x6^3 + x6^4,
        x5^5 + x5^4*x6 + x5^3*x6^2 + x5^2*x6^3 + x5*x6^4 + x6^5,
        x6^6 + 2147483646
    ]
    #! format: on

    ku = Groebner.Examples.ku10(k=GF(2^31 - 1), internal_ordering=:degrevlex)
    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = gens(parent(first(ku)))
    gb = Groebner.groebner(ku, reduced=true)

    @test gb == [
        x9 + 1272065637 * x10 + 875418006,
        x8 + 1529540685 * x10 + 617942964,
        x7 + 1539832471 * x10 + 607651173,
        x6 + 1314302432 * x10 + 833181218,
        x5 + 1453635454 * x10 + 693848197,
        x4 + 673118236 * x10 + 1474365406,
        x3 + 269783061 * x10 + 1877700587,
        x2 + 1042807874 * x10 + 1104675778,
        x1 + 389079675 * x10 + 1758403970,
        x10^2 + 1222705397 * x10 + 924778249
    ]
end

@testset "groebner certify" begin
    R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"], internal_ordering=:lex)

    @test Groebner.groebner([x, y], certify=true) == [y, x]
    @test Groebner.groebner([y, x], certify=true) == [y, x]

    fs = [x^2 + y, x * y]
    @test Groebner.groebner(fs, certify=true) == [y^2, x * y, x^2 + y]

    root = Groebner.Examples.rootn(3, k=QQ, internal_ordering=:degrevlex)
    x1, x2, x3 = gens(parent(first(root)))
    gb = Groebner.groebner(root, certify=true)
    @test gb == [x1 + x2 + x3, x2^2 + x2 * x3 + x3^2, x3^3 - 1]

    root = Groebner.Examples.rootn(6, k=QQ, internal_ordering=:deglex)
    x1, x2, x3, x4, x5, x6 = gens(parent(first(root)))
    gb = Groebner.groebner(root, certify=true)
    #! format: off
    @test gb == [
        x1 + x2 + x3 + x4 + x5 + x6,
        x2^2 + x2*x3 + x3^2 + x2*x4 + x3*x4 + x4^2 + x2*x5 + x3*x5 + x4*x5 + x5^2 + x2*x6 + x3*x6 + x4*x6 + x5*x6 + x6^2,
        x3^3 + x3^2*x4 + x3*x4^2 + x4^3 + x3^2*x5 + x3*x4*x5 + x4^2*x5 + x3*x5^2 + x4*x5^2 + x5^3 + x3^2*x6 + x3*x4*x6 + x4^2*x6 + x3*x5*x6 + x4*x5*x6 + x5^2*x6 + x3*x6^2 + x4*x6^2 + x5*x6^2 + x6^3,
        x4^4 + x4^3*x5 + x4^2*x5^2 + x4*x5^3 + x5^4 + x4^3*x6 + x4^2*x5*x6 + x4*x5^2*x6 + x5^3*x6 + x4^2*x6^2 + x4*x5*x6^2 + x5^2*x6^2 + x4*x6^3 + x5*x6^3 + x6^4,
        x5^5 + x5^4*x6 + x5^3*x6^2 + x5^2*x6^3 + x5*x6^4 + x6^5,
        x6^6 - 1
    ]
    #! format: on

    ku = Groebner.Examples.ku10(k=QQ, internal_ordering=:degrevlex)
    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = gens(parent(first(ku)))
    gb = Groebner.groebner(ku, certify=true)

    @test gb == [
        x9 + 230965 // 116706 * x10 - 697789 // 116706,
        x8 + 92386 // 335511 * x10 + 578636 // 335511,
        x7 + 5173616 // 8563059 * x10 - 30862793 // 8563059,
        x6 + 87951472 // 192015891 * x10 + 488096201 // 192015891,
        x5 + 131980 // 18759 * x10 - 56944 // 18759,
        x4 + 32995 // 6084 * x10 - 63415 // 6084,
        x3 - 263960 // 229671 * x10 + 493631 // 229671,
        x2 - 65990 // 17043 * x10 + 151205 // 17043,
        x1 - 26396 // 3273 * x10 + 19850 // 3273,
        x10^2 - 71551 // 26396 * x10 + 45155 // 26396
    ]

    system = [
        x1 + BigInt(10)^100 // BigInt(7)^50 * x2,
        BigInt(90) * x1 * x2 + BigInt(19)^50 + x10^3,
        BigInt(2^31 - 1) * x1^2 + BigInt(2^30 + 2)^10 * x2^2
    ]
    @test Groebner.groebner(system) == Groebner.groebner(system, certify=true)
end

@testset "groebner orderings" begin
    R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], internal_ordering=:deglex)

    gb = Groebner.groebner([x, y, z, w], ordering=Groebner.Lex(y, x, w, z))
    @test gb == [z, w, x, y]

    # Parent ring persists
    for aa_ord in [:lex, :deglex, :degrevlex]
        R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], internal_ordering=aa_ord)

        @test_throws DomainError Groebner.groebner([x, y, z, w], ordering=Groebner.Lex(x))
        @test_throws DomainError Groebner.groebner([x, y], ordering=Groebner.DegLex(x, y))
        @test_throws DomainError Groebner.groebner([x, y], ordering=Groebner.DegRevLex(x, y))

        for gb_ord in [
            Groebner.Lex(),
            Groebner.DegLex(),
            Groebner.DegRevLex(),
            Groebner.Lex(y, w, z, x),
            Groebner.DegLex(y, w, z, x),
            Groebner.DegRevLex(y, w, z, x),
            Groebner.WeightedOrdering(Dict([x, y, z, w] .=> [1, 1, 1, 1])),
            Groebner.Lex(x, y) * Groebner.DegLex(z, w),
            Groebner.Lex(x) * Groebner.DegLex(y) * Groebner.DegRevLex(z, w),
            Groebner.MatrixOrdering([x, y, z, w], [1 2 3 4; 5 6 7 8])
        ]
            gb0 = Groebner.groebner([R(0), R(0)], ordering=gb_ord)
            gb1 = Groebner.groebner([R(5)], ordering=gb_ord)
            gb2 = Groebner.groebner([x + y + z], ordering=gb_ord)
            @test R == parent(gb0[1]) == parent(gb1[1]) == parent(gb2[1])
            @test gb0 == [R(0)] && gb1 == [R(1)] && gb2 == [x + y + z]
        end
    end

    # Normalization is correct
    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
    f = [7x + 3, 11y + 5z + 7]
    @test Groebner.groebner(f) == [y + (5 // 11) * z + 7 // 11, x + 3 // 7]
    @test Groebner.groebner(f, ordering=Groebner.Lex(y, x, z)) ==
          [x + 3 // 7, y + (5 // 11) * z + 7 // 11]
    @test Groebner.groebner(f, ordering=Groebner.Lex(z, x, y)) ==
          [x + 3 // 7, z + (11 // 5) * y + 7 // 5]

    f = [7x + 3y + 5z + 11]
    @test Groebner.groebner(f) == [x + (3 // 7) * y + (5 // 7) * z + 11 // 7]
    @test Groebner.groebner(f, ordering=Groebner.Lex(y, x, z)) ==
          [(7 // 3) * x + y + (5 // 3) * z + 11 // 3]
    @test Groebner.groebner(f, ordering=Groebner.Lex(z, x, y)) ==
          [(7 // 5) * x + (3 // 5) * y + z + 11 // 5]

    # The order of output polynomials is correct
    R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], internal_ordering=:deglex)
    for i in 1:10
        xs = Random.shuffle([x, y, z, w])
        polys = [x + 1, y + 2, z + 3, w - 4]

        gb = Groebner.groebner(Random.shuffle(polys), ordering=Groebner.Lex(xs))
        @test map(leading_term, gb) == reverse(xs)

        gb = Groebner.groebner(Random.shuffle(polys), ordering=Groebner.DegLex(xs))
        @test map(leading_term, gb) == reverse(xs)

        gb = Groebner.groebner(Random.shuffle(polys), ordering=Groebner.DegRevLex(xs))
        @test map(leading_term, gb) == reverse(xs)
    end

    # Handles rings that are not the same
    R, (a, b) = polynomial_ring(QQ, ["a", "b"])
    R, (A, B, C) = polynomial_ring(QQ, ["a", "b", "c"])
    @test groebner([a, b], ordering=Lex(B, A)) == [a, b]
    @test groebner([a, b], ordering=Lex(B, A, C)) == [a, b]
    @test groebner([a, b], ordering=DegLex(B) * DegRevLex(A)) == [a, b]
    @test groebner([a, b], ordering=DegLex(A) * DegRevLex(B)) == [b, a]

    # Correctness of Lex, DegLex, DegRevLex comparing with AbstractAlgebra
    for aa_ord in [:lex, :deglex, :degrevlex]
        R, (x1, x2, x3, x4) =
            polynomial_ring(QQ, ["x1", "x2", "x3", "x4"], internal_ordering=aa_ord)
        cases = [
            [x1, x2, x3, x4],
            [x2 * x3 + 2x4, x1 * x2 + 3x3],
            [x1 + 2x2 + 3x3 + 4x4, x1, x2, (x1 + x2 + x3 + x4)^3, x3, x4],
            Groebner.Examples.rootn(5)
        ]
        for case in cases
            xs = gens(parent(case[1]))
            for vars in Combinatorics.permutations(xs)
                for (gb_ord, _quot_aa_ord) in [
                    (Groebner.Lex(vars), Meta.quot(:lex)),
                    (Groebner.DegLex(vars), Meta.quot(:deglex)),
                    (Groebner.DegRevLex(vars), Meta.quot(:degrevlex))
                ]
                    gb = Groebner.groebner(Random.shuffle(case), ordering=gb_ord)

                    # A hack to get the same ranking of variables as in Groebner
                    eval(
                        :(
                            (__R, _) = polynomial_ring(
                                QQ,
                                map(repr, $(vars)),
                                internal_ordering=$_quot_aa_ord
                            )
                        )
                    )
                    _xs = sort(gens(__R), by=repr)

                    _case = map(poly -> evaluate(poly, _xs), case)
                    _gb = Groebner.groebner(_case)
                    gb = map(poly -> evaluate(poly, _xs), gb)

                    @test gb == _gb
                end
            end
        end
    end

    # WeightedOrdering
    R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], internal_ordering=:deglex)
    @test_throws DomainError Groebner.groebner(
        [x, y],
        ordering=Groebner.WeightedOrdering(Dict([x, y, z] .=> [1, 0, 2]))
    )
    @test_throws DomainError Groebner.groebner(
        [x, y],
        ordering=Groebner.WeightedOrdering(Dict([x, y, z, x, x] .=> [1, 0, 1, 9, 10]))
    )

    R, (x1, x2, x3, x4, x5, x6) = QQ["x1", "x2", "x3", "x4", "x5", "x6"]
    @test_throws DomainError Groebner.groebner(
        [x1],
        ordering=Groebner.WeightedOrdering(Dict([x1, x1, x1, x1, x1, x1] .=> [-1, 0, 0, 0, 0, 0]))
    )

    # ProductOrdering
    ord = Groebner.Lex(x6, x2) * Groebner.Lex(x4, x1, x3)
    @test Groebner.groebner([x1, x2], ordering=ord) == [x1, x2]
    @test Groebner.groebner([x4, x6, x3], ordering=ord) == [x3, x4, x6]

    ord = Groebner.Lex(x6, x2, x5) * Groebner.Lex(x4, x1, x3)
    @test [x3, x1, x4, x5, x2, x6] == Groebner.groebner([x1, x2, x3, x4, x5, x6], ordering=ord)
    f1, f2 = x4 + x1 + x3, x1^2
    F = [f1 + x6 * x2 * x5, x6 * x2 * x5, f1 - f2]
    @test Groebner.groebner(F, ordering=ord)[1:2] == [f2, f1]

    ord = Groebner.Lex(x6, x2, x1, x5) * Groebner.Lex(x4, x1, x3)
    @test [x3, x4, x5, x1, x2, x6] == Groebner.groebner([x1, x2, x3, x4, x5, x6], ordering=ord)

    # MatrixOrdering
    c = Groebner.Examples.cyclicn(4)
    x = gens(parent(c[1]))
    ord = Groebner.MatrixOrdering(
        x,
        [
            1 0 0 0
            0 1 0 0
            0 0 1 0
            0 0 0 1
        ]
    )
    @test Groebner.groebner(c, ordering=ord) == Groebner.groebner(c, ordering=Groebner.Lex())

    ord = Groebner.MatrixOrdering(
        x,
        [
            0 0 0 1
            0 0 1 0
            0 1 0 0
            1 0 0 0
        ]
    )
    @test Groebner.groebner(c, ordering=ord) ==
          Groebner.groebner(c, ordering=Groebner.Lex(reverse(x)))

    ord = Groebner.MatrixOrdering(
        x,
        [
            1 1 1 1
            1 0 0 0
            0 1 0 0
            0 0 1 0
            0 0 0 1
        ]
    )
    @test Groebner.groebner(c, ordering=ord) == Groebner.groebner(c, ordering=Groebner.DegLex())

    ord = Groebner.MatrixOrdering(
        x,
        [
            1 1 1 1
            0 0 0 -1
            0 0 -1 0
            0 -1 0 0
            -1 0 0 0
        ]
    )
    @test Groebner.groebner(c, ordering=ord) == Groebner.groebner(c, ordering=Groebner.DegRevLex())
end

@testset "groebner parent rings" begin
    R, x = polynomial_ring(QQ, "x")
    @test R == parent(first(Groebner.groebner([x], ordering=Groebner.Lex())))
    @test R == parent(first(Groebner.groebner([x], ordering=Groebner.DegLex())))
    @test R == parent(first(Groebner.groebner([x], ordering=Groebner.DegRevLex())))

    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], internal_ordering=:deglex)

    fs = [y^2, x]
    gb = Groebner.groebner(fs, ordering=Groebner.Lex())
    @test parent(gb[1]) == R
    @test Groebner.groebner(fs, ordering=Groebner.Lex()) == [y^2, x]

    gb = Groebner.groebner(fs, ordering=Groebner.DegRevLex())
    @test parent(gb[1]) == R
    @test gb == [x, y^2]

    gb = Groebner.groebner(fs, ordering=Groebner.DegLex())
    @test parent(gb[1]) == R
    @test gb == [x, y^2]

    noon = Groebner.Examples.noonn(2, internal_ordering=:lex)
    gb = Groebner.groebner(noon, ordering=Groebner.Lex())
    x1, x2 = gens(parent(first(noon)))
    @test gb == [
        x2^5 - 10 // 11 * x2^4 - 11 // 5 * x2^3 + 2 * x2^2 + 331 // 1100 * x2 - 11 // 10,
        x1 + x2^4 - 10 // 11 * x2^3 - 11 // 10 * x2^2 + x2 - 10 // 11
    ]

    gb = Groebner.groebner(noon, ordering=Groebner.DegRevLex())
    x1, x2 = gens(parent(first(noon)))
    @test gb == [
        x1^2 - x2^2 - 10 // 11 * x1 + 10 // 11 * x2,
        x2^3 + 10 // 11 * x1 * x2 - 10 // 11 * x2^2 - 11 // 10 * x2 + 1,
        x1 * x2^2 - 11 // 10 * x1 + 1
    ]

    gb = Groebner.groebner(noon, ordering=Groebner.DegLex())
    x1, x2 = gens(parent(first(noon)))
    @test gb == [
        x1^2 - x2^2 - 10 // 11 * x1 + 10 // 11 * x2,
        x2^3 + 10 // 11 * x1 * x2 - 10 // 11 * x2^2 - 11 // 10 * x2 + 1,
        x1 * x2^2 - 11 // 10 * x1 + 1
    ]
end

@testset "groebner monoms" begin
    for domain in (GF(2^31 - 1), QQ)
        for system in [
            Groebner.Examples.cyclicn(2, k=domain),
            Groebner.Examples.noonn(4, k=domain, internal_ordering=:degrevlex),
            Groebner.Examples.katsuran(5, k=domain, internal_ordering=:degrevlex),
            Groebner.Examples.kinema(k=domain, internal_ordering=:degrevlex)
        ]
            results = []
            for monoms in [:dense, :packed]
                gb = Groebner.groebner(system, monoms=monoms)
                push!(results, gb)
                @test Groebner.isgroebner(gb)
            end
            @test length(unique(results)) == 1
        end
    end
end

@testset "groebner zeros" begin
    for field in [GF(2), GF(2^31 - 1), QQ]
        R, (x, y, z) = polynomial_ring(field, ["x", "y", "z"], internal_ordering=:lex)

        @test Groebner.groebner([x - x]) == [R(0)]
        @test Groebner.groebner([R(0), R(0), R(0)]) == Groebner.groebner([R(0)]) == [R(0)]
        @test Groebner.groebner([x, R(0), y, R(0)]) == [y, x]

        @test_throws DomainError Groebner.groebner([])
    end
end

@testset "isgroebner zeros" begin
    for field in [GF(2), GF(2^31 - 1), QQ]
        R, (x, y, z) = polynomial_ring(field, ["x", "y", "z"], internal_ordering=:lex)

        @test Groebner.isgroebner([x - x])
        @test Groebner.isgroebner([R(0), R(0), R(0)])
        @test Groebner.isgroebner([R(0), R(0), R(1)])
        @test Groebner.isgroebner([x, R(0), y, R(0)])

        @test_throws DomainError Groebner.isgroebner([])
    end
end

@testset "normalform zeros" begin
    for field in [GF(2), GF(2^31 - 1), QQ]
        R, (x, y, z) = polynomial_ring(field, ["x", "y", "z"], internal_ordering=:lex)

        @test Groebner.normalform([R(0)], [R(0)]) == [R(0)]
        @test Groebner.normalform([R(0), R(0), R(0)], [x + 2, x, x + 1]) == [x + 2, x, x + 1]
        @test Groebner.normalform([R(0), R(0), R(0)], R(0)) == R(0)

        @test_throws DomainError Groebner.normalform([], x)
        @test_throws DomainError Groebner.normalform([], [x])
    end
end

@testset "normalform checks" begin
    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], internal_ordering=:lex)

    for check in [false, true]
        @test Groebner.normalform([x, y], y, check=check) == R(0)
    end
    @test_throws DomainError Groebner.normalform([x, x + 1], y, check=true)
    @test_throws DomainError Groebner.normalform([x, x + 1], [y], check=true)
end

@testset "groebner arithmetic" begin
    R, (x, y, z) = polynomial_ring(GF(10007), ["x", "y", "z"], internal_ordering=:degrevlex)

    for arithmetic in [:auto, :delayed, :signed, :basic]
        @test Groebner.groebner([x, y], arithmetic=arithmetic, linalg=:deterministic) ==
              Groebner.groebner([y, x]) ==
              [y, x]

        sys = Groebner.Examples.rootn(3, k=GF(10007), internal_ordering=:degrevlex)
        x1, x2, x3 = gens(parent(first(sys)))
        gb1 = Groebner.groebner(sys, arithmetic=arithmetic, linalg=:deterministic)
        gb2 = Groebner.groebner(sys)
        @test gb1 == gb2 == [x1 + x2 + x3, x2^2 + x2 * x3 + x3^2, x3^3 + 10006]

        sys = Groebner.Examples.cyclicn(7, k=GF(2^25 - 39))
        gb1 = Groebner.groebner(sys, arithmetic=arithmetic, linalg=:deterministic)
        gb2 = Groebner.groebner(sys)
        @test gb1 == gb2

        sys = Groebner.Examples.noonn(7, k=GF(2^25 - 39))
        gb1 = Groebner.groebner(sys, arithmetic=arithmetic, linalg=:deterministic)
        gb2 = Groebner.groebner(sys)
        @test gb1 == gb2
    end
end

@testset "groebner linear algebra" begin
    R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"], internal_ordering=:lex)

    for linalg in [:deterministic, :randomized]
        @test Groebner.groebner([x, y], linalg=linalg) == Groebner.groebner([y, x]) == [y, x]

        fs = [x^2 + y, x * y]
        @test Groebner.groebner(fs, linalg=linalg) == Groebner.groebner(fs) == [y^2, x * y, x^2 + y]

        root = Groebner.Examples.rootn(3, k=GF(2^31 - 1), internal_ordering=:degrevlex)
        x1, x2, x3 = gens(parent(first(root)))
        gb1 = Groebner.groebner(root, linalg=linalg)
        gb2 = Groebner.groebner(root)
        @test gb1 == gb2 == [x1 + x2 + x3, x2^2 + x2 * x3 + x3^2, x3^3 + 2147483646]

        root = Groebner.Examples.rootn(6, k=GF(2^31 - 1), internal_ordering=:degrevlex)
        x1, x2, x3, x4, x5, x6 = gens(parent(first(root)))
        gb1 = Groebner.groebner(root, linalg=linalg)
        gb2 = Groebner.groebner(root)
        @test gb1 == gb2

        ku = Groebner.Examples.ku10(k=GF(2^31 - 1), internal_ordering=:degrevlex)
        x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = gens(parent(first(ku)))
        gb1 = Groebner.groebner(ku)
        gb2 = Groebner.groebner(ku, linalg=linalg)

        @test gb1 ==
              gb2 ==
              [
                  x9 + 1272065637 * x10 + 875418006,
                  x8 + 1529540685 * x10 + 617942964,
                  x7 + 1539832471 * x10 + 607651173,
                  x6 + 1314302432 * x10 + 833181218,
                  x5 + 1453635454 * x10 + 693848197,
                  x4 + 673118236 * x10 + 1474365406,
                  x3 + 269783061 * x10 + 1877700587,
                  x2 + 1042807874 * x10 + 1104675778,
                  x1 + 389079675 * x10 + 1758403970,
                  x10^2 + 1222705397 * x10 + 924778249
              ]
    end
end

@testset "groebner modular-hard problems" begin
    function get_test_system1(R, N)
        (x1, x2, x3, x4) = AbstractAlgebra.gens(R)
        system = [
            x1 + x2 + x3 + x4,
            x1 * x2 + x1 * x3 + x1 * x4 + x2 * x3 + x2 * x4 + x3 * x4,
            x1 * x2 * x3 + x1 * x2 * x4 + x1 * x3 * x4 + x2 * x3 * x4,
            x1 * x2 * x3 * x4 + N
        ]
        result = [
            x1 + x2 + x3 + x4,
            x2^2 + x2 * x3 + x3^2 + x2 * x4 + x3 * x4 + x4^2,
            x3^3 + x3^2 * x4 + x3 * x4^2 + x4^3,
            x4^4 - N
        ]
        system, result
    end

    R, (x1, x2, x3, x4) =
        polynomial_ring(QQ, ["x1", "x2", "x3", "x4"], internal_ordering=:degrevlex)

    N = prod(map(BigInt, nextprimes(2^30 + 3, 5)))
    system, result = get_test_system1(R, N)
    # this should take about 10 primes
    gb = Groebner.groebner(system)
    @test gb == result

    N = prod(map(BigInt, nextprimes(2^31 - 1, 5)))
    system, result = get_test_system1(R, N)
    # this should take about 10 primes
    gb = Groebner.groebner(system)
    @test gb == result

    N = prod(map(BigInt, nextprimes(2^31 - 1, 100)))
    system, result = get_test_system1(R, N)
    # around 200 primes are required
    gb = Groebner.groebner(system)
    @test gb == result

    N = prod(map(BigInt, nextprimes(2^31 - 1, 5_000)))
    # around 10k primes are required
    system, result = get_test_system1(R, N)
    gb = Groebner.groebner(system)
    @test gb == result

    system = [x1 - (2^31 - 1) * x2 - (2^30 + 3) * x3]
    @test Groebner.groebner(system) == system

    for start_of_range in [2^10, 2^20, 2^30, 2^40]
        for size_of_range in [10, 100, 1000]
            N = prod(Primes.nextprimes(BigInt(start_of_range), size_of_range))
            system = [x1 + N // (N + 1), x3 - (N + 1) // (N - 1), x2 + 42]
            gb = Groebner.groebner(system)
            @test gb == [x3 - (N + 1) // (N - 1), x2 + 42, x1 + N // (N + 1)]
        end
    end
end

@testset "groebner strange example" begin
    get_rand_poly(x, d, n) = sum([rand(1:5) * prod(x .^ rand(0:d, length(x))) for _ in 1:n])

    n = 10
    R, x = polynomial_ring(GF(17), ["x$i" for i in 1:n], internal_ordering=:degrevlex)
    f1 = get_rand_poly(x, 10, 2^15)
    @test Groebner.groebner([f1]) == [divexact(f1, leading_coefficient(f1))]

    n = 6
    R, x = polynomial_ring(GF(17), ["x$i" for i in 1:n], internal_ordering=:degrevlex)
    f1 = get_rand_poly(x, 3, 2^10)
    f2 = get_rand_poly(x, 2, 2^10)
    gb = Groebner.groebner([f1, f2])
    @info "GB contains polynomials of lengths: $(sort(map(length, gb)))"
end

@testset "groebner many variables" begin
    function n_variable_set(n)
        R, x = polynomial_ring(QQ, ["x$i" for i in 1:n])
        f = [sum(prod(x[i:(n - k)], init=1) for i in 1:(k + 1)) for k in 0:(n - 1)]
        f
    end

    function test_n_variables(n)
        f = n_variable_set(n)
        x = gens(parent(f[1]))
        gb = Groebner.groebner(f)

        evencf(i) = isone(i) ? 0 // 1 : (2(i - 1)) // (2i - 1)
        oddcf(i) = isone(i) ? 0 // 1 : (2(i - 1) - 1) // (2(i - 1))
        cf(i, n) = (iseven(n) ? oddcf(i) : evencf(i)) * (-1n)^(i == div(n, 2) + 1)

        ans = [x[div(n, 2) - i + 2] - cf(i, n) for i in 1:(div(n, 2) + 1)]

        @test Groebner.groebner(f, ordering=Groebner.DegRevLex()) == ans
        @test Groebner.isgroebner(gb)
        @test all(iszero, Groebner.normalform(gb, f))
        @test gb == ans
    end

    # up to 63
    test_n_variables(8)
    test_n_variables(13)
    test_n_variables(16)
    test_n_variables(20)
    test_n_variables(32)
    for n in 32:8:63
        test_n_variables(n)
    end

    # up to 127
    for n in [64, 100, 101, 127]
        @info "Variables:" n
        test_n_variables(n)
    end

    # up to 257
    for n in [128, 256, 257]
        @info "Variables:" n
        R, x = polynomial_ring(QQ, ["x$i" for i in 1:n])
        f = x
        Groebner.groebner(f)
    end
end

# For Lex, DegLex, and DegRevLex, we can have total degrees up to 2^31
@testset "groebner large exponents" begin
    R, (x, y) = polynomial_ring(QQ, ["x", "y"], internal_ordering=:degrevlex)

    # up to 2^8-1
    for (i, d) in enumerate(4:2:255)
        f = [x^d - 1, x * y + 2]
        m, n, k = div(d, 2), div(d, 2) + 1, div(d, 2) - 1
        gb = Groebner.groebner(f)
        @test gb == [
            x * y + 2,
            x^m + (-1)^(i) // BigInt(2)^m * y^m,
            y^n + (-1)^(i + 1) * BigInt(2)^n * x^k
        ]
    end

    # up to 5^6 < 2^14
    for i in 1:6
        u, v = 3^i, 5^i
        f = [x^u * y^v - 1, x^v + y^u]
        gb = Groebner.groebner(f)
        @test gb == [x^v + y^u, y^(v + u) + x^(v - u), x^u * y^v - 1]
    end

    # above 2^16-1
    f = [x^(2^16) + y]
    @test f == Groebner.groebner(f)

    # up to 5^13 < 2^32
    for i in 7:13
        u, v = 3^i, 5^i
        f = [x^u * y^v - 1, x^v + y^u]
        gb = Groebner.groebner(f)
        @test gb == [x^v + y^u, y^(v + u) + x^(v - u), x^u * y^v - 1]
    end

    # total degree 2^30
    f = [x^1073741824 + y]
    @test f == Groebner.groebner(f)

    # total degree 2^31
    f = [x^1073741824 * y^1073741824 + y]
    @test f == Groebner.groebner(f)

    # this should fail.
    # f = [x^4294967295 * y^4294967295 + y]
    # @test f == Groebner.groebner(f)

    for i in [8, 16, 32, 64]
        u, v = i >> 1, i
        f = [x^u * y^v + y^8, x^8 - y^8]
        gb = Groebner.groebner(f)
    end
end

@testset "homogenization, basic" begin
    # TODO: broken for char. 113
    for field in [GF(2^10 + 7), GF(2^62 + 135), QQ]
        R, x = polynomial_ring(field, "x")
        @test Groebner.groebner([x^2 - 1, (x + 1) * x^2], homogenize=:yes) == [x + 1]
        @test Groebner.groebner([x^2 - 1, (x + 1) * x^2], homogenize=:no) == [x + 1]
        @test Groebner.groebner([x^2 - 1, (x + 1) * x^2], homogenize=:auto) == [x + 1]

        for ordering in [:lex, :deglex, :degrevlex]
            R, (x, y) = polynomial_ring(field, ["x", "y"], internal_ordering=ordering)

            @test Groebner.groebner([R(3)], homogenize=:yes) == [R(1)]
            @test Groebner.groebner([R(0)], homogenize=:yes) == [R(0)]
            @test Groebner.groebner([x], homogenize=:yes) == [x]
            @test Groebner.groebner([x + 11y], homogenize=:yes) == [x + 11y]
            @test Groebner.groebner([x, y], homogenize=:yes) == [y, x]
            @test Groebner.isgroebner(Groebner.groebner([x + 1, y + 2], homogenize=:yes))

            if ordering === :lex
                gb = Groebner.groebner([x + y^2, x^2 - 1], homogenize=:yes)
                @test (y^4 - 1) in gb && (x + y^2) in gb
            end

            for case in [
                Groebner.Examples.katsuran(2, k=field, internal_ordering=ordering),
                Groebner.Examples.katsuran(3, k=field, internal_ordering=ordering),
                Groebner.Examples.katsuran(4, k=field, internal_ordering=ordering),
                Groebner.Examples.noonn(2, k=field, internal_ordering=ordering),
                Groebner.Examples.noonn(3, k=field, internal_ordering=ordering)
            ]
                gb = Groebner.groebner(case, homogenize=:yes)
                @test Groebner.isgroebner(gb)
            end
        end
    end
end

@testset "homogenization, orderings" begin
    for field in [GF(113), GF(2^62 + 135), QQ]
        # Test that the basis obtained with the use of homogenization
        # *coincides* with the one obtained without it
        R, (x, y, z) = polynomial_ring(field, ["x", "y", "z"])
        for case in [
            (system=[R(1)], ord=Groebner.Lex()),
            (system=[R(5)], ord=Groebner.DegRevLex()),
            (system=[x, y, z], ord=Groebner.Lex()),
            (system=[x, y, z], ord=Groebner.Lex(z) * Groebner.Lex(x, y)),
            (system=[x^2 + 1, x * y + 2, y * z + 3], ord=Groebner.DegLex()),
            (
                system=[x^20 - x^15 - x^5 + y * z, x * y^10 + x^3 * y^3 + x * y],
                ord=Groebner.DegRevLex()
            ),
            (system=[x + 5y, x + 7z, y + 11z], ord=Groebner.Lex(z) * Groebner.Lex(x, y)),
            (system=[x + 5y + 11z], ord=Groebner.DegRevLex(y) * Groebner.Lex(z, x)),
            (system=[x + 5y + 11z], ord=Groebner.DegRevLex(z) * Groebner.DegRevLex(y, x)),
            (system=[x * y^2 + x + 1, y * z^2 + y + 1, z^4 - z^2 - 1], ord=Groebner.Lex()),
            (system=[x * y^2 + x + 1, y * z^2 + y + 1, z^4 - z^2 - 1], ord=Groebner.DegRevLex()),
            (
                system=[x * y^2 + x + 1, y * z^2 + y + 1, z^4 - z^2 - 1],
                ord=Groebner.DegRevLex(z) * Groebner.DegRevLex(x, y)
            ),
            (
                system=[x * y^2 + x + 1, y * z^2 + y + 1, x * y * z^4 - x * z^2 - y + z],
                ord=Groebner.DegRevLex(y) * Groebner.DegRevLex(z, x)
            ),
            (system=Groebner.Examples.katsuran(4, k=field), ord=Groebner.DegRevLex()),
            (system=Groebner.Examples.rootn(5, k=field), ord=Groebner.Lex()),
            (system=Groebner.Examples.reimern(5, k=field), ord=Groebner.DegRevLex())
        ]
            ord = case.ord
            system = case.system
            @debug "" system ord field

            gb1 = Groebner.groebner(system, ordering=ord, homogenize=:no)
            gb2 = Groebner.groebner(system, ordering=ord, homogenize=:yes)

            @test Groebner.isgroebner(gb1, ordering=ord)

            @test gb1 == gb2

            # Also test learn / apply
            field == QQ && continue

            context3, gb3 = Groebner.groebner_learn(system, ordering=ord, homogenize=:yes)
            context4, gb4 = Groebner.groebner_learn(system, ordering=ord, homogenize=:no)
            for _ in 1:4
                flag3, gb33 = Groebner.groebner_apply!(context3, system)
                flag4, gb44 = Groebner.groebner_apply!(context4, system)
                @test flag3 && flag4
                @test gb1 == gb3 == gb4 == gb33 == gb44
            end
        end
    end
end

@testset "groebner, change matrix" begin
    R, (x, y, z) = polynomial_ring(GF(2^30 + 3), ["x", "y", "z"], internal_ordering=:degrevlex)
    f = [x * y * z - 1, x * y + x * z + y * z, x + y + z]
    g, m = Groebner.groebner_with_change_matrix(f)
    @test m * f == g
    @test size(m) == (length(g), length(f))
    @test Groebner.isgroebner(g)

    for ground in [GF(2^30 + 3), QQ]
        R, (x, y) = polynomial_ring(ground, ["x", "y"], internal_ordering=:degrevlex)

        cases = [
            [R(0), R(0), R(0)],
            [R(0), R(1), R(0)],
            [x],
            [y],
            [x, x, x, x, R(0), R(0), x, x, x, x],
            [x + y, R(1), x^5 - y^5],
            [x^200 + y^250, x^100 * y^200 + 1],
            Groebner.Examples.katsuran(4, k=ground),
            Groebner.Examples.noonn(4, k=ground),
            Groebner.Examples.cyclicn(4, k=ground),
            [x^i * y + x for i in 1:500]
        ]

        for f in cases
            g, m = Groebner.groebner_with_change_matrix(f)
            @test m * f == g
            @test size(m) == (length(g), length(f))
            @test Groebner.isgroebner(g)
        end
    end

    f = Groebner.Examples.katsuran(4, k=QQ, internal_ordering=:lex)
    g, m = Groebner.groebner_with_change_matrix(f, ordering=Groebner.DegRevLex())
    @test m * f == g
    @test_throws DomainError Groebner.groebner_with_change_matrix(f, ordering=Groebner.DegLex())
    @test_throws DomainError Groebner.groebner_with_change_matrix(f)
end

@testset "groebner, multi-threading, Zp" begin
    @info "Testing multi-threading over Zp using $(nthreads()) threads"

    R, (x, y, z) = GF(2^31 - 1)["x", "y", "z"]

    s = [R(1), R(2), R(4)]
    @test Groebner.groebner(s, threaded=:yes) ==
          Groebner.groebner(s, threaded=:no) ==
          Groebner.groebner(s, threaded=:auto) ==
          [R(1)]

    s = [x^(2^10) + 1, y^3000 + 2, z^1000 + 3]
    @test Groebner.groebner(s, threaded=:yes) ==
          Groebner.groebner(s, threaded=:no) ==
          Groebner.groebner(s, threaded=:auto) ==
          [z^1000 + 3, y^3000 + 2, x^(2^10) + 1]

    linalg = [:deterministic, :randomized]
    threaded = [:yes, :no, :auto]
    grid = [(linalg=l, threaded=t) for l in linalg for t in threaded]
    for system in [
        Groebner.Examples.katsuran(3, internal_ordering=:lex, k=GF(2^31 - 1)),
        Groebner.Examples.katsuran(4, internal_ordering=:lex, k=GF(2^20 + 7)),
        Groebner.Examples.katsuran(7, internal_ordering=:degrevlex, k=GF(2^31 - 1)),
        Groebner.Examples.eco5(internal_ordering=:deglex, k=GF(2^20 + 7)),
        Groebner.Examples.eco5(internal_ordering=:degrevlex, k=GF(2^20 + 7)),
        Groebner.Examples.cyclicn(5, internal_ordering=:degrevlex, k=GF(2^20 + 7)),
        Groebner.Examples.cyclicn(6, internal_ordering=:degrevlex, k=GF(2^20 + 7)),
        Groebner.Examples.ojika4(internal_ordering=:lex, k=GF(2^20 + 7)),
        Groebner.Examples.henrion5(internal_ordering=:degrevlex, k=GF(2^20 + 7)),
        Groebner.Examples.noonn(6, internal_ordering=:degrevlex, k=GF(2^10 + 7)),
        Groebner.Examples.ku10(internal_ordering=:degrevlex, k=GF(2^10 + 7))
    ]
        results = Array{Any}(undef, size(grid))
        for (i, kw) in enumerate(grid)
            gb = Groebner.groebner(system; kw...)
            results[i] = gb
        end
        @test length(unique(results)) == 1
        @test Groebner.isgroebner(results[1])
    end
end

@testset "groebner, multi-threading, QQ" begin
    @info "Testing multi-threading over QQ using $(nthreads()) threads"

    R, (x, y) = QQ["x", "y"]

    threaded = [:yes, :no, :auto]
    grid = [(threaded=t,) for t in threaded]
    for system in [
        [x - 1, y - 2],
        [x + (BigInt(2)^1000 + 1) // 2^61, x * y + BigInt(2)^(2^10)],
        Groebner.Examples.katsuran(3, internal_ordering=:lex, k=QQ),
        Groebner.Examples.katsuran(4, internal_ordering=:lex, k=QQ),
        Groebner.Examples.eco5(internal_ordering=:degrevlex, k=QQ),
        Groebner.Examples.econ(4, internal_ordering=:degrevlex, k=QQ),
        Groebner.Examples.cyclicn(5, internal_ordering=:degrevlex, k=QQ),
        Groebner.Examples.ojika4(internal_ordering=:lex, k=QQ),
        Groebner.Examples.noonn(6, internal_ordering=:degrevlex, k=QQ),
        Groebner.Examples.ku10(internal_ordering=:degrevlex, k=QQ),
        [x + BigInt(2)^1_000 // (BigInt(2)^1_000 - 1), y - BigInt(2)^1_000],
        [x + BigInt(2)^10_000 // (BigInt(2)^10_000 - 1), y^2 - BigInt(2)^10_000 * y]
    ]
        results = Array{Any}(undef, size(grid))
        for (i, kw) in enumerate(grid)
            gb = Groebner.groebner(system; kw...)
            results[i] = gb
        end
        @test length(unique(results)) == 1
        @test Groebner.isgroebner(results[1])
    end
end

function test_params(
    rng,
    nvariables,
    maxdegs,
    nterms,
    npolys,
    grounds,
    coeffssize,
    orderings,
    linalgs,
    monoms,
    homogenizes
)
    boot = 1
    for _ in 1:boot
        for nv in nvariables
            for md in maxdegs
                for nt in nterms
                    for np in npolys
                        for k in grounds
                            for ord in orderings
                                for cf in coeffssize
                                    if ord == :lex
                                        md = min(2, md)
                                    end
                                    set = Groebner.Examples.random_generating_set(
                                        rng,
                                        k,
                                        ord,
                                        nv,
                                        md,
                                        nt,
                                        np,
                                        cf
                                    )
                                    isempty(set) && continue

                                    for linalg in linalgs
                                        for monom in monoms
                                            for homogenize in homogenizes
                                                try
                                                    gb = Groebner.groebner(
                                                        set,
                                                        linalg=linalg,
                                                        monoms=monom,
                                                        homogenize=homogenize
                                                    )
                                                    flag = Groebner.isgroebner(gb)
                                                    if !flag
                                                        @error "Beda!" nv md nt np k ord monom
                                                        println("Rng:\n", rng)
                                                        println("Set:\n", set)
                                                        println("Gb:\n", gb)
                                                    end
                                                    @test flag
                                                catch err
                                                    @error "Beda!" nv md nt np k ord monom
                                                    println(err)
                                                    println("Rng:\n", rng)
                                                    println("Set:\n", set)
                                                    rethrow(err)
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

@testset "groebner random stress tests" begin
    rng = Random.MersenneTwister(42)

    nvariables = [2, 3]
    maxdegs    = [2, 4]
    nterms     = [1, 2, 4]
    npolys     = [1, 4, 100]
    grounds    = [GF(1031), GF(2^50 + 55), AbstractAlgebra.QQ]
    coeffssize = [3, 1000, 2^31 - 1, BigInt(2)^80]
    orderings  = [:degrevlex, :lex, :deglex]
    linalgs    = [:deterministic, :randomized]
    monoms     = [:auto, :dense, :packed]
    homogenize = [:yes, :auto]
    p          = prod(map(length, (nvariables, maxdegs, nterms, npolys, grounds, orderings, coeffssize, linalgs, monoms, homogenize)))
    @info "Producing $p small random tests for groebner. This may take a minute"

    test_params(
        rng,
        nvariables,
        maxdegs,
        nterms,
        npolys,
        grounds,
        coeffssize,
        orderings,
        linalgs,
        monoms,
        homogenize
    )
end
