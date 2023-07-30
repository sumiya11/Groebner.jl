using AbstractAlgebra

@testset "groebner simple" begin
    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:degrevlex)

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

    R, (x1, x2) = PolynomialRing(GF(2^31 - 1), ["x1", "x2"], ordering=:degrevlex)

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

@testset "groebner noreduce" begin
    R, (x, y) = PolynomialRing(GF(2^31 - 1), ["x", "y"], ordering=:lex)

    fs = [x + y^2, x * y - y^2]
    gb = Groebner.groebner(fs, reduced=false)
    @test gb ≂ [x + y^2, x * y + 2147483646 * y^2, y^3 + y^2]

    fs = [x + y, x^2 + y]
    gb = Groebner.groebner(fs, reduced=false)
    @test gb ≂ [x + y, x^2 + y, y^2 + y]

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
    R, (x,) = PolynomialRing(QQ, ["x"], ordering=:degrevlex)
    @test Groebner.groebner([x]) == [x]
    @test Groebner.groebner([4x]) == [x]
    @test Groebner.groebner([x + 1]) == [x + 1]
    @test Groebner.groebner([4x + 1]) == [x + 1 // 4]
    @test Groebner.groebner([x + 5]) == [x + 5]
    @test Groebner.groebner([x + 2^10]) == [x + 2^10]
    @test Groebner.groebner([x + 2^30]) == [x + 2^30]
    @test Groebner.groebner([x + 2^30 + 3]) == [x + 2^30 + 3]
    @test Groebner.groebner([x + 2^31 - 1]) == [x + 2^31 - 1]
    @test Groebner.groebner([x + 2^31 - 1, x^2]) == [1]
    @test Groebner.groebner([(3232323 // 7777)x + 7777 // 3232323]) ==
          [x + 60481729 // 10447911976329]
    @test Groebner.groebner([((2^31 - 1) // 1)x + 1]) == [x + 1 // 2147483647]
    @test Groebner.groebner([(1 // (2^31 - 1))x + 1]) == [x + 2147483647]
    @test Groebner.groebner([1 // (2^30 + 3) * x^2 + (2^30 + 3)x + 1 // (1073741831)]) ==
          [x^2 + 1152921511049297929x + (2^30 + 3) // 1073741831]

    R, (x, y) = PolynomialRing(QQ, ["x", "y"], ordering=:degrevlex)

    fs = [x, y]
    G = Groebner.groebner(fs)
    @test G ≂ [y, x]

    fs = [5y, 3 // 4 * x]
    G = Groebner.groebner(fs)
    @test G ≂ [y, x]

    fs = [x^2 - (1 // 6) * y, x * y]
    G = Groebner.groebner(fs)
    @test G ≂ [y^2, x * y, x^2 - (1 // 6)y]

    fs = [
        QQ(11, 3) * x^2 * y - QQ(2, 4) * x * y - QQ(1, 7) * y,
        QQ(1, 1) * x * y + QQ(7, 13) * x
    ]
    G = Groebner.groebner(fs)
    @test G ≂ [y^2 + 7 // 13 * y, x * y + 7 // 13 * x, x^2 - 3 // 22 * x + 39 // 539 * y]

    root = Groebner.rootn(3, ground=QQ, ordering=:degrevlex)
    gb = Groebner.groebner(root)
    @test Groebner.isgroebner(gb)

    root = Groebner.rootn(4, ground=QQ, ordering=:degrevlex)
    gb = Groebner.groebner(root)
    @test Groebner.isgroebner(gb)

    root = Groebner.rootn(8, ground=QQ, ordering=:degrevlex)
    gb = Groebner.groebner(root)
    @test Groebner.isgroebner(gb)

    noon = Groebner.noonn(3, ground=QQ, ordering=:degrevlex)
    gb = Groebner.groebner(noon)
    @test Groebner.isgroebner(gb)

    noon = Groebner.noonn(4, ground=QQ, ordering=:degrevlex)
    gb = Groebner.groebner(noon)
    @test Groebner.isgroebner(gb)
end

@testset "groebner corners cases" begin
    R, (x, y, z, w) = PolynomialRing(QQ, ["x", "y", "z", "w"], ordering=:degrevlex)

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

    R, (x, y, z, w) = PolynomialRing(QQ, ["x", "y", "z", "w"], ordering=:lex)
    fs = [(12345678 // 12347)x, (222222221111123 // 2131232232097)y + z]
    G = Groebner.groebner(fs)
    @test G == [y + 2131232232097 // 222222221111123 * z, x]
end

@testset "groebner output sorted" begin
    for K in [GF(11), GF(307), GF(12007), GF(2^31 - 1), QQ]
        R, (x, y, z) = PolynomialRing(K, ["x", "y", "z"], ordering=:lex)

        @test Groebner.groebner([x, y]) == Groebner.groebner([y, x]) == [y, x]

        # TODO: find out why this fails on 1.6
        @static if VERSION >= v"1.7.0"
            @test Groebner.groebner([x^2 + y, x * y]) == [y^2, x * y, x^2 + y]
        end

        @test Groebner.groebner([3x + 2, 5y]) == [y, x + K(2) // K(3)]
    end
end

@testset "monomial overflow" begin
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"], ordering=:degrevlex)
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

@testset "groebner reduced basis" begin
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
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"], ordering=:lex)

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
    R, (x, y, z, w) = PolynomialRing(QQ, ["x", "y", "z", "w"], ordering=:deglex)

    Groebner.groebner([x, y], ordering=Groebner.Lex())
    Groebner.groebner([x, y], ordering=Groebner.DegLex())
    Groebner.groebner([x, y], ordering=Groebner.DegRevLex())

    # Test that the order of output polynomials is correct
    for i in 1:10
        xs = shuffle([x, y, z, w])
        polys = [x + 1, y + 2, z + 3, w - 4]

        gb = Groebner.groebner(shuffle(polys), ordering=Groebner.Lex(xs))
        @test map(string ∘ leading_term, gb) == map(string, reverse(xs))

        gb = Groebner.groebner(shuffle(polys), ordering=Groebner.DegLex(xs))
        @test map(string ∘ leading_term, gb) == map(string, reverse(xs))

        gb = Groebner.groebner(shuffle(polys), ordering=Groebner.DegRevLex(xs))
        @test map(string ∘ leading_term, gb) == map(string, reverse(xs))
    end

    # for aa_ord in [:lex, :deglex, :degrevlex]
    #     R, (x1, x2, x3, x4) = PolynomialRing(QQ, ["x1", "x2", "x3", "x4"], ordering=aa_ord)
    #     cases = []
    # end

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
    ord = Groebner.Lex(x6, x2, x5) * Groebner.Lex(x4, x1, x3)
    @test [x3, x1, x4, x5, x2, x6] ==
          Groebner.groebner([x1, x2, x3, x4, x5, x6], ordering=ord)

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
    @test_throws AssertionError Groebner.groebner([x1, x2], ordering=ord)
end

@testset "groebner on-the-fly order change" begin
    R, x = PolynomialRing(QQ, "x")
    @test R == parent(Groebner.groebner([x], ordering=Groebner.Lex())[1])
    @test R == parent(Groebner.groebner([x], ordering=Groebner.DegLex())[1])
    @test R == parent(Groebner.groebner([x], ordering=Groebner.DegRevLex())[1])

    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"], ordering=:deglex)

    fs = [y^2, x]
    x, y, z = gens(parent(first(Groebner.groebner(fs, ordering=Groebner.Lex()))))
    @test AbstractAlgebra.ordering(parent(x)) == :lex
    @test Groebner.groebner(fs, ordering=Groebner.Lex()) == [y^2, x]
    @test AbstractAlgebra.ordering(
        parent(first(Groebner.groebner(fs, ordering=Groebner.Lex())))
    ) === :lex

    x, y, z = gens(parent(first(Groebner.groebner(fs, ordering=Groebner.DegRevLex()))))
    @test Groebner.groebner(fs, ordering=Groebner.DegRevLex()) == [x, y^2]
    @test AbstractAlgebra.ordering(
        parent(first(Groebner.groebner(fs, ordering=Groebner.DegRevLex())))
    ) === :degrevlex

    x, y, z = gens(parent(first(fs)))
    @test Groebner.groebner(fs, ordering=Groebner.DegLex()) == [x, y^2]
    @test AbstractAlgebra.ordering(
        parent(first(Groebner.groebner(fs, ordering=Groebner.DegLex())))
    ) === :deglex

    noon = Groebner.change_ordering(Groebner.noonn(2), :lex)
    gb = Groebner.groebner(noon, ordering=Groebner.Lex())
    x1, x2 = gens(parent(first(gb)))
    @test gb == [
        x2^5 - 10 // 11 * x2^4 - 11 // 5 * x2^3 + 2 * x2^2 + 331 // 1100 * x2 - 11 // 10,
        x1 + x2^4 - 10 // 11 * x2^3 - 11 // 10 * x2^2 + x2 - 10 // 11
    ]

    gb = Groebner.groebner(noon, ordering=Groebner.DegRevLex())
    x1, x2 = gens(parent(first(gb)))
    @test gb == [
        x1^2 - x2^2 - 10 // 11 * x1 + 10 // 11 * x2,
        x2^3 + 10 // 11 * x1 * x2 - 10 // 11 * x2^2 - 11 // 10 * x2 + 1,
        x1 * x2^2 - 11 // 10 * x1 + 1
    ]

    gb = Groebner.groebner(noon, ordering=Groebner.DegLex())
    x1, x2 = gens(parent(first(gb)))
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
        R, (x, y, z) = PolynomialRing(field, ["x", "y", "z"], ordering=:lex)

        @test Groebner.groebner([x - x]) == [R(0)]
        @test Groebner.groebner([R(0), R(0), R(0)]) == Groebner.groebner([R(0)]) == [R(0)]
        @test Groebner.groebner([x, R(0), y, R(0)]) == [y, x]

        @test_throws AssertionError Groebner.groebner([])
    end
end

@testset "isgroebner zeros" begin
    for field in [GF(2), GF(2^31 - 1), QQ]
        R, (x, y, z) = PolynomialRing(field, ["x", "y", "z"], ordering=:lex)

        @test Groebner.isgroebner([x - x])
        @test Groebner.isgroebner([R(0), R(0), R(0)])
        @test Groebner.isgroebner([R(0), R(0), R(1)])
        @test Groebner.isgroebner([x, R(0), y, R(0)])

        @test_throws AssertionError Groebner.isgroebner([])
    end
end

@testset "normalform zeros" begin
    for field in [GF(2), GF(2^31 - 1), QQ]
        R, (x, y, z) = PolynomialRing(field, ["x", "y", "z"], ordering=:lex)

        @test Groebner.normalform([R(0)], [R(0)]) == [R(0)]
        @test Groebner.normalform([R(0), R(0), R(0)], [x + 2, x, x + 1]) ==
              [x + 2, x, x + 1]
        @test Groebner.normalform([R(0), R(0), R(0)], R(0)) == R(0)

        @test_throws AssertionError Groebner.normalform([], x)
        @test_throws AssertionError Groebner.normalform([], [x])
    end
end

@testset "normalform checks" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"], ordering=:lex)

    for check in [false, true]
        @test Groebner.normalform([x, y], y, check=check) == R(0)
    end
    @test_throws DomainError Groebner.normalform([x, x + 1], y, check=true)
    @test_throws DomainError Groebner.normalform([x, x + 1], [y], check=true)
end

@testset "kbase checks" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"], ordering=:lex)

    for check in [false, true]
        b = Groebner.kbase([x, y, z], check=check)
        @test b == [R(1)]
    end
    @test_throws DomainError Groebner.kbase([x, x + 1], check=true)
    @test_throws DomainError Groebner.kbase([x, x + 1], check=true)
end

@testset "groebner linear algebra" begin
    R, (x, y, z) = PolynomialRing(GF(2^31 - 1), ["x", "y", "z"], ordering=:lex)

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
