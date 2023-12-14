using AbstractAlgebra
using Combinatorics
import Random

@testset "groebner simple" begin
    R, (x, y) = polynomial_ring(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)

    @test_throws DomainError Groebner.groebner([])
    @test_throws DomainError Groebner.groebner(Vector{elem_type(R)}())

    fs = [x, x, x]
    @test Groebner.groebner(fs) ≂ [x]

    fs = [x, x, y, y, y, x, x, y]
    @test Groebner.groebner(fs) ≂ [x, y]

    fs = [R(1)]
    @test Groebner.groebner(fs) ≂ [R(1)]

    fs = [x^2 + y, y * x - 1, R(1), y^4]
    @test Groebner.groebner(fs) ≂ [R(1)]

    fs = [2x, 3y, 4x + 5y]
    @test Groebner.groebner(fs) ≂ [x, y]

    R, (x1, x2) = polynomial_ring(GF(2^31 - 1), ["x1", "x2"], ordering=:degrevlex)

    fs = [
        1 * x1^2 * x2^2 + 2 * x1^2 * x2,
        1 * x1^2 * x2^2 + 3 * x1^2 * x2 + 5 * x1 * x2^2,
        1 * x1^2 * x2^2
    ]
    gb = Groebner.groebner(fs)
    @test gb ≂ [x1 * x2^2, x1^2 * x2]

    fs = [x1 * x2^2 + x1 * x2, x1^2 * x2^2 + 2 * x1 * x2, 2 * x1^2 * x2^2 + x1 * x2]
    gb = Groebner.groebner(fs)
    @test gb ≂ [x1 * x2]
end

@testset "groebner no autoreduce" begin
    R, (x, y) = polynomial_ring(GF(2^31 - 1), ["x", "y"], ordering=:lex)

    fs = [x + y^2, x * y - y^2]
    gb = Groebner.groebner(fs, reduced=false)
    @test_broken gb ≂ [x + y^2, x * y + 2147483646 * y^2, y^3 + y^2]
    @test Groebner.isgroebner(gb)

    fs = [x + y, x^2 + y]
    gb = Groebner.groebner(fs, reduced=false)
    @test_broken gb ≂ [x + y, x^2 + y, y^2 + y]
    @test Groebner.isgroebner(gb)

    fs = [y, x * y + x]
    gb = Groebner.groebner(fs, reduced=false)
    @test gb ≂ [x, y]

    fs = [x^2 + 5, 2y^2 + 3]
    gb = Groebner.groebner(fs, reduced=false)
    @test gb ≂ [y^2 + 1073741825, x^2 + 5]

    fs = [y^2 + x, x^2 * y + y, y + 1]
    gb = Groebner.groebner(fs, reduced=false)
    @test gb == [R(1)]

    fs = [x^2 * y^2, 2 * x * y^2 + 3 * x * y]
    gb = Groebner.groebner(fs, reduced=false)
    @test gb ≂ [x * y^2 + 1073741825 * x * y, x^2 * y]
end

@testset "groebner modular" begin
    for modular in [:classic_modular, :learn_and_apply]
        R, (x,) = polynomial_ring(QQ, ["x"], ordering=:degrevlex)
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
        @test Groebner.groebner([((2^31 - 1) // 1)x + 1], modular=modular) ==
              [x + 1 // 2147483647]
        @test Groebner.groebner([(1 // (2^31 - 1))x + 1], modular=modular) ==
              [x + 2147483647]
        @test Groebner.groebner(
            [1 // (2^30 + 3) * x^2 + (2^30 + 3)x + 1 // (1073741831)],
            modular=modular
        ) == [x^2 + 1152921511049297929x + (2^30 + 3) // 1073741831]

        R, (x, y) = polynomial_ring(QQ, ["x", "y"], ordering=:degrevlex)

        fs = [x, y]
        G = Groebner.groebner(fs, modular=modular)
        @test G ≂ [y, x]

        fs = [5y, 3 // 4 * x]
        G = Groebner.groebner(fs, modular=modular)
        @test G ≂ [y, x]

        fs = [x^2 - (1 // 6) * y, x * y]
        G = Groebner.groebner(fs, modular=modular)
        @test G ≂ [y^2, x * y, x^2 - (1 // 6)y]

        fs = [
            QQ(11, 3) * x^2 * y - QQ(2, 4) * x * y - QQ(1, 7) * y,
            QQ(1, 1) * x * y + QQ(7, 13) * x
        ]
        G = Groebner.groebner(fs, modular=modular)
        @test G ≂
              [y^2 + 7 // 13 * y, x * y + 7 // 13 * x, x^2 - 3 // 22 * x + 39 // 539 * y]

        root = Groebner.rootn(3, ground=QQ, ordering=:degrevlex)
        gb = Groebner.groebner(root, modular=modular)
        @test Groebner.isgroebner(gb)

        root = Groebner.rootn(4, ground=QQ, ordering=:degrevlex)
        gb = Groebner.groebner(root, modular=modular)
        @test Groebner.isgroebner(gb)

        root = Groebner.rootn(8, ground=QQ, ordering=:degrevlex)
        gb = Groebner.groebner(root, modular=modular)
        @test Groebner.isgroebner(gb)

        noon = Groebner.noonn(3, ground=QQ, ordering=:degrevlex)
        gb = Groebner.groebner(noon, modular=modular)
        @test Groebner.isgroebner(gb)

        noon = Groebner.noonn(4, ground=QQ, ordering=:degrevlex)
        gb = Groebner.groebner(noon, modular=modular)
        @test Groebner.isgroebner(gb)

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
end

@testset "groebner corner cases" begin
    R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], ordering=:degrevlex)

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

    @test G ≂ [
        w^3,
        z * w + 1692834523553 // 4063974 * w^2,
        z^2 - 16935085031076 // 16933225 * w^2,
        y - 12345 // 4106996 * z - 4115 // 1026749 * w,
        x + 6172839 // 2053498 * z + 4115226 // 1026749 * w
    ]

    G = [(2^30 + 3)x, (1 // (2^30 + 3))y]
    G = Groebner.groebner(G)
    @test G == [y, x]

    G = [(2^31 - 1) * x + y, (1 // (2^31 - 1)) * y + 1]
    G = Groebner.groebner(G)
    @test G == [y + 2147483647, x - 1]

    R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], ordering=:lex)
    fs = [(12345678 // 12347)x, (222222221111123 // 2131232232097)y + z]
    G = Groebner.groebner(fs)
    @test G == [y + 2131232232097 // 222222221111123 * z, x]
end

@testset "groebner output sorted" begin
    for K in [GF(11), GF(307), GF(12007), GF(2^31 - 1), QQ]
        R, (x, y, z) = polynomial_ring(K, ["x", "y", "z"], ordering=:lex)

        @test Groebner.groebner([x, y]) == Groebner.groebner([y, x]) == [y, x]

        # TODO: find out why this fails on 1.6
        @static if VERSION >= v"1.7.0"
            @test Groebner.groebner([x^2 + y, x * y]) == [y^2, x * y, x^2 + y]
        end

        @test Groebner.groebner([3x + 2, 5y]) == [y, x + K(2) // K(3)]
    end
end

@testset "monomial overflow" begin
    R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"], ordering=:degrevlex)
    for monoms in [:auto, :dense, :packed, :sparse]
        gb_1 = [x * y^100 + y, x^100 * y + y^100, y^199 + 2147483646 * x^99 * y]
        gb_2 = [x * y^200 + y, x^200 * y + y^200, y^399 + 2147483646 * x^199 * y]
        gb_3 = [x * y^1000 + y, x^1000 * y + y^1000, y^1999 + 2147483646 * x^999 * y]
        @test Groebner.groebner([x^100 * y + y^100, x * y^100 + y], monoms=monoms) == gb_1
        @test Groebner.groebner([x^200 * y + y^200, x * y^200 + y], monoms=monoms) == gb_2
        @test Groebner.groebner([x^1000 * y + y^1000, x * y^1000 + y], monoms=monoms) ==
              gb_3

        @test Groebner.isgroebner(gb_1)
        @test Groebner.isgroebner(gb_2)
        @test Groebner.isgroebner(gb_3)

        @test Groebner.normalform(gb_1, [x, y, R(1), R(0), x^1000]) ==
              [x, y, R(1), R(0), x^1000]
        @test Groebner.normalform(gb_2, [x, y, R(1), R(0), x^1000]) ==
              [x, y, R(1), R(0), x^1000]
        @test Groebner.normalform(gb_3, [x, y, R(1), R(0), x^10]) ==
              [x, y, R(1), R(0), x^10]
        @test Groebner.normalform(gb_3, [x, y, R(1), R(0), x^1000]) ==
              [x, y, R(1), R(0), x^1000]
    end
end

@testset "groebner autoreduce" begin
    root = Groebner.rootn(3, ground=GF(2^31 - 1), ordering=:degrevlex)
    x1, x2, x3 = gens(parent(first(root)))
    gb = Groebner.groebner(root, reduced=true)
    @test gb == [x1 + x2 + x3, x2^2 + x2 * x3 + x3^2, x3^3 + 2147483646]

    root = Groebner.rootn(6, ground=GF(2^31 - 1), ordering=:degrevlex)
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

    ku = Groebner.ku10(ground=GF(2^31 - 1), ordering=:degrevlex)
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
    R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"], ordering=:lex)

    @test Groebner.groebner([x, y], certify=true) == [y, x]
    @test Groebner.groebner([y, x], certify=true) == [y, x]

    fs = [x^2 + y, x * y]
    @test Groebner.groebner(fs, certify=true) == [y^2, x * y, x^2 + y]

    root = Groebner.rootn(3, ground=QQ, ordering=:degrevlex)
    x1, x2, x3 = gens(parent(first(root)))
    gb = Groebner.groebner(root, certify=true)
    @test gb == [x1 + x2 + x3, x2^2 + x2 * x3 + x3^2, x3^3 - 1]

    root = Groebner.rootn(6, ground=QQ, ordering=:deglex)
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

    ku = Groebner.ku10(ground=QQ, ordering=:degrevlex)
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

@testset "groebner maxpairs" begin
    # TODO: fix some bugs in maxpairs
    s = Groebner.noonn(5, ordering=:degrevlex, ground=GF(2^31 - 1))
    gb = Groebner.groebner(s)
    # @test gb == Groebner.groebner(s, maxpairs=100)
    # @test gb == Groebner.groebner(s, maxpairs=10)
    # @test gb == Groebner.groebner(s, maxpairs=2)
    # @test gb == Groebner.groebner(s, maxpairs=1)
    @test_throws AssertionError Groebner.groebner(s, maxpairs=0)

    s = Groebner.katsuran(5, ordering=:degrevlex, ground=QQ)
    gb = Groebner.groebner(s)
    # @test gb == Groebner.groebner(s, maxpairs=100)
    # @test gb == Groebner.groebner(s, maxpairs=10)
    # @test gb == Groebner.groebner(s, maxpairs=2)
    # @test gb == Groebner.groebner(s, maxpairs=1)
    @test_throws AssertionError Groebner.groebner(s, maxpairs=0)
end

@testset "groebner orderings" begin
    R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], ordering=:deglex)

    gb = Groebner.groebner([x, y, z, w], ordering=Groebner.Lex(y, x, w, z))
    @test gb == [z, w, x, y]

    # Parent ring persists
    for aa_ord in [:lex, :deglex, :degrevlex]
        R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], ordering=aa_ord)

        @test_throws DomainError Groebner.groebner([x, y, z, w], ordering=Groebner.Lex(x))
        @test_throws DomainError Groebner.groebner([x, y], ordering=Groebner.DegLex(x, y))
        @test_throws DomainError Groebner.groebner(
            [x, y],
            ordering=Groebner.DegRevLex(x, y)
        )

        for gb_ord in [
            Groebner.Lex(),
            Groebner.DegLex(),
            Groebner.DegRevLex(),
            Groebner.Lex(y, w, z, x),
            Groebner.DegRevLex(y, w, z, x),
            Groebner.DegRevLex(y, w, z, x),
            Groebner.WeightedOrdering([1, 1, 1, 1]),
            Groebner.Lex(x, y) * Groebner.DegLex(z, w),
            Groebner.Lex(x) * Groebner.DegLex(y) * Groebner.DegRevLex(z, w),
            Groebner.MatrixOrdering([1 2 3 4; 5 6 7 8])
        ]
            gb0 = Groebner.groebner([R(0), R(0)], ordering=gb_ord)
            gb1 = Groebner.groebner([R(5)], ordering=gb_ord)
            gb2 = Groebner.groebner([x + y + z], ordering=gb_ord)
            @test R == parent(gb0[1]) == parent(gb1[1]) == parent(gb2[1])
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
    R, (x, y, z, w) = polynomial_ring(QQ, ["x", "y", "z", "w"], ordering=:deglex)
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

    # Correctness of Lex, DegLex, DegRevLex comparing with AbstractAlgebra
    for aa_ord in [:lex, :deglex, :degrevlex]
        R, (x1, x2, x3, x4) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4"], ordering=aa_ord)
        cases = [
            [x1, x2, x3, x4],
            [x2 * x3 + 2x4, x1 * x2 + 3x3],
            [x1 + 2x2 + 3x3 + 4x4, x1, x2, (x1 + x2 + x3 + x4)^3, x3, x4],
            Groebner.rootn(5)
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
                                ordering=$_quot_aa_ord
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
    @test_throws DomainError Groebner.groebner(
        [x, y],
        ordering=Groebner.WeightedOrdering([1, 0])
    )
    @test_throws DomainError Groebner.groebner(
        [x, y],
        ordering=Groebner.WeightedOrdering([1, 0, 1, 9, 10])
    )

    R, (x1, x2, x3, x4, x5, x6) = QQ["x1", "x2", "x3", "x4", "x5", "x6"]
    @test_throws AssertionError Groebner.groebner(
        [x1],
        ordering=Groebner.WeightedOrdering([-1, 0, 0, 0, 0, 0])
    )

    # ProductOrdering
    ord = Groebner.Lex(x6, x2) * Groebner.Lex(x4, x1, x3)
    @test_throws DomainError Groebner.groebner([x], ordering=ord)

    ord = Groebner.Lex(x6, x2, x5) * Groebner.Lex(x4, x1, x3)
    @test [x3, x1, x4, x5, x2, x6] ==
          Groebner.groebner([x1, x2, x3, x4, x5, x6], ordering=ord)

    ord = Groebner.Lex(x6, x2, x1, x5) * Groebner.Lex(x4, x1, x3)
    @test [x3, x4, x5, x1, x2, x6] ==
          Groebner.groebner([x1, x2, x3, x4, x5, x6], ordering=ord)

    # MatrixOrdering
    ord = Groebner.MatrixOrdering([
        1 0 0 0 1 2
        0 1 0 0 -2 -1
        0 0 1 0 0 0
        0 0 0 1 0 0
    ])
    Groebner.groebner([x1, x2], ordering=ord)

    ord = Groebner.MatrixOrdering([
        1 0 0 0
        0 1 0 0
        0 0 1 0
    ])
    @test_throws DomainError Groebner.groebner([x1, x2], ordering=ord)
end

@testset "groebner parent rings" begin
    R, x = polynomial_ring(QQ, "x")
    @test R == parent(first(Groebner.groebner([x], ordering=Groebner.Lex())))
    @test R == parent(first(Groebner.groebner([x], ordering=Groebner.DegLex())))
    @test R == parent(first(Groebner.groebner([x], ordering=Groebner.DegRevLex())))

    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], ordering=:deglex)

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

    noon = Groebner.noonn(2, ordering=:lex)
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
            Groebner.cyclicn(2, ground=domain),
            Groebner.noonn(4, ground=domain, ordering=:degrevlex),
            Groebner.katsuran(5, ground=domain, ordering=:degrevlex),
            Groebner.kinema(ground=domain, ordering=:degrevlex)
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
        R, (x, y, z) = polynomial_ring(field, ["x", "y", "z"], ordering=:lex)

        @test Groebner.groebner([x - x]) == [R(0)]
        @test Groebner.groebner([R(0), R(0), R(0)]) == Groebner.groebner([R(0)]) == [R(0)]
        @test Groebner.groebner([x, R(0), y, R(0)]) == [y, x]

        @test_throws DomainError Groebner.groebner([])
    end
end

@testset "isgroebner zeros" begin
    for field in [GF(2), GF(2^31 - 1), QQ]
        R, (x, y, z) = polynomial_ring(field, ["x", "y", "z"], ordering=:lex)

        @test Groebner.isgroebner([x - x])
        @test Groebner.isgroebner([R(0), R(0), R(0)])
        @test Groebner.isgroebner([R(0), R(0), R(1)])
        @test Groebner.isgroebner([x, R(0), y, R(0)])

        @test_throws DomainError Groebner.isgroebner([])
    end
end

@testset "normalform zeros" begin
    for field in [GF(2), GF(2^31 - 1), QQ]
        R, (x, y, z) = polynomial_ring(field, ["x", "y", "z"], ordering=:lex)

        @test Groebner.normalform([R(0)], [R(0)]) == [R(0)]
        @test Groebner.normalform([R(0), R(0), R(0)], [x + 2, x, x + 1]) ==
              [x + 2, x, x + 1]
        @test Groebner.normalform([R(0), R(0), R(0)], R(0)) == R(0)

        @test_throws DomainError Groebner.normalform([], x)
        @test_throws DomainError Groebner.normalform([], [x])
    end
end

@testset "normalform checks" begin
    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], ordering=:lex)

    for check in [false, true]
        @test Groebner.normalform([x, y], y, check=check) == R(0)
    end
    @test_throws DomainError Groebner.normalform([x, x + 1], y, check=true)
    @test_throws DomainError Groebner.normalform([x, x + 1], [y], check=true)
end

@testset "kbase checks" begin
    R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"], ordering=:lex)

    for check in [false, true]
        b = Groebner.kbase([x, y, z], check=check)
        @test b == [R(1)]
    end
    @test_throws DomainError Groebner.kbase([x, x + 1], check=true)
    @test_throws DomainError Groebner.kbase([x, x + 1], check=true)
end

@testset "groebner linear algebra" begin
    R, (x, y, z) = polynomial_ring(GF(2^31 - 1), ["x", "y", "z"], ordering=:lex)

    for linalg in [:deterministic, :randomized]
        @test Groebner.groebner([x, y], linalg=linalg) ==
              Groebner.groebner([y, x]) ==
              [y, x]

        fs = [x^2 + y, x * y]
        @test Groebner.groebner(fs, linalg=linalg) ==
              Groebner.groebner(fs) ==
              [y^2, x * y, x^2 + y]

        root = Groebner.rootn(3, ground=GF(2^31 - 1), ordering=:degrevlex)
        x1, x2, x3 = gens(parent(first(root)))
        gb1 = Groebner.groebner(root, linalg=linalg)
        gb2 = Groebner.groebner(root)
        @test gb1 == gb2 == [x1 + x2 + x3, x2^2 + x2 * x3 + x3^2, x3^3 + 2147483646]

        root = Groebner.rootn(6, ground=GF(2^31 - 1), ordering=:degrevlex)
        x1, x2, x3, x4, x5, x6 = gens(parent(first(root)))
        gb1 = Groebner.groebner(root, linalg=linalg)
        gb2 = Groebner.groebner(root)
        @test gb1 == gb2

        ku = Groebner.ku10(ground=GF(2^31 - 1), ordering=:degrevlex)
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
