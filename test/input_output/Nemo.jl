using Test, Nemo, Groebner

@testset "Nemo.jl, univariate" begin
    for ff in [Nemo.Native.GF(2^62 + 135), Nemo.GF(2^62 + 135)]
        R, x = Nemo.polynomial_ring(ff, "x")
        @test Groebner.groebner([R(2)]) == [R(1)]
        @test Groebner.groebner([R(0), R(0)]) == [R(0)]
        @test Groebner.groebner([R(0), R(3), R(0)]) == [R(1)]
        @test Groebner.groebner([R(0), R(3), R(0)]) == [R(1)]
        @test Groebner.groebner([x + 100]) == [x + 100]
    end

    for ff in [
        Nemo.Native.GF(2^31 - 1),
        Nemo.Native.GF(2^62 + 135),
        Nemo.GF(2^62 + 135),
        Nemo.GF(Nemo.ZZRingElem(2^62 + 135)),
        Nemo.QQ
    ]
        for gb_ord in [Groebner.Lex(), Groebner.DegLex(), Groebner.DegRevLex()]
            R, x = Nemo.polynomial_ring(ff, "x")
            @test Groebner.groebner([x^2 - 4, x + 2], ordering=gb_ord) == [x + 2]
        end
    end
end

@testset "Nemo.jl, generic" begin
    # Test generic coefficients with some interesting fields from Nemo.jl.

    # Single extension
    K, a = Nemo.finite_field(3, 2, "a")
    R, (X, Y) = K["X", "Y"]
    @test groebner([a * X - Y, X * Y - 1]) == [Y^2 + 2 * a, X + (2 * a + 1) * Y]

    # TODO: Tower of extensions
    # K, a = Nemo.finite_field(3, 2, "a")
    # Kx, x = Nemo.polynomial_ring(K, "x")
    # L, b = Nemo.finite_field(x^3 + x^2 + x + 2, "b")
    # R, (X, Y) = L["X", "Y"]
    # sys = [a * b * X - Y, X * Y - 1]
    # res = [Y^2 + b^5 + b^4 + 2*b^2 + 2*b + 2, X + (b^3 + 2*b^2 + b)*Y]
    # @test groebner(sys) == res

    # Q bar
    Q_bar = Nemo.algebraic_closure(Nemo.QQ)
    R, (X, Y) = Q_bar["X", "Y"]
    e2 = Nemo.root_of_unity(Q_bar, 5, 2)
    e3 = Nemo.root_of_unity(Q_bar, 5, 3)
    e4 = Nemo.root_of_unity(Q_bar, 5, 4)
    @test groebner([e3 * X - e2 * Y]) == [X - e4 * Y]

    # Cyclic extension
    K, a = Nemo.cyclotomic_field(5)
    R, (X, Y) = K["X", "Y"]
    @test groebner([a * X + 1]) == [X - a^3 - a^2 - a - 1]
    @test isgroebner([X - a^3 - a^2 - a - 1])
    @test normalform([X - a^3 - a^2 - a - 1], [X, X - a^3 - a^2 - a - 1]) ==
          [R(a^3 + a^2 + a + 1), R(0)]
    @test_throws DomainError groebner_with_change_matrix(
        [X - a^3 - a^2 - a - 1],
        ordering=Groebner.DegRevLex()
    )
end

@testset "Nemo.jl, input-output" begin
    # R, (x, y) = Nemo.GF(Nemo.ZZRingElem(Primes.nextprime(BigInt(2)^100)))["x", "y"]
    # @test_throws DomainError Groebner.groebner([x, y])

    # R, (x, y) = Nemo.GF(2, 2)["x", "y"]
    # @test_throws DomainError Groebner.groebner([x, y])

    R, (x, y) = Nemo.ZZ["x", "y"]
    @test_throws DomainError Groebner.groebner([x, y])

    R_, (a) = Nemo.QQ["a"]
    R, (x, y) = R_["x", "y"]
    @test_throws DomainError Groebner.groebner([x, y])

    # R_, (a) = Nemo.QQ["a"]
    # R, (x, y) = AbstractAlgebra.FractionField(R_)["x", "y"]
    # @test_throws DomainError Groebner.groebner([x, y])

    nemo_orderings_to_test = [:lex, :deglex, :degrevlex]
    nemo_grounds_to_test = [
        Nemo.Native.GF(2^62 + 135),
        Nemo.Native.GF(2^31 - 1),
        Nemo.Native.GF(17),
        Nemo.GF(2^62 + 135),
        Nemo.GF(2^31 - 1),
        Nemo.QQ
    ]

    for ord in nemo_orderings_to_test
        for ground in nemo_grounds_to_test
            R, x = Nemo.polynomial_ring(ground, "x")
            gb = Groebner.groebner([(x - 1) * (x + 8), (x + 8) * (x + 10)])
            @test gb == [(x + 8)]

            R, (x, y) = Nemo.polynomial_ring(ground, ["x", "y"], internal_ordering=ord)
            fs = [x^2 * y + 3, (2^31 - 5) * x - (2^31 - 4) * y]
            gb = Groebner.groebner(fs)
            @test parent(gb[1]) == R
            @test Groebner.isgroebner(gb)
        end
    end

    c = Groebner.Examples.cyclicn(6, k=Nemo.GF(2^30 + 3))
    gb1 = Groebner.groebner(c)
    trace, gb2 = Groebner.groebner_learn(c)
    flag, gb3 = Groebner.groebner_apply!(trace, c)
    flag, _gb4 = Groebner.groebner_apply!(trace, (c, c, c, c))
    @test gb1 == gb2 == gb3 == _gb4[1]

    c = Groebner.Examples.cyclicn(6, k=Nemo.Native.GF(2^30 + 3))
    gb1 = Groebner.groebner(c)
    trace, gb2 = Groebner.groebner_learn(c)
    flag, gb3 = Groebner.groebner_apply!(trace, c)
    flag, _gb4 = Groebner.groebner_apply!(trace, (c, c, c, c))
    @test gb1 == gb2 == gb3 == _gb4[1]

    # Test for different Groebner.jl orderings
    for nemo_ord in [:lex, :deglex, :degrevlex]
        for ground in nemo_grounds_to_test
            R, (x,) = Nemo.polynomial_ring(ground, ["x"], internal_ordering=nemo_ord)
            for gb_ord in [
                Groebner.Lex(),
                Groebner.DegLex(),
                Groebner.DegRevLex(),
                Groebner.Lex(x),
                Groebner.DegLex(x)
            ]
                gb = Groebner.groebner([x^2], ordering=gb_ord)
                @test parent(first(gb)) == R
                @test gb == [x^2]
            end

            R, (x, y) = Nemo.polynomial_ring(ground, ["x", "y"], internal_ordering=nemo_ord)
            fs = [x^2 + 3, y - 1]
            for gb_ord in [
                Groebner.Lex(),
                Groebner.DegLex(),
                Groebner.DegRevLex(),
                Groebner.Lex(x, y),
                Groebner.Lex(y, x)
            ]
                gb = Groebner.groebner(fs, ordering=gb_ord)
                @test parent(first(gb)) == R
                @test all(in(fs), gb) && all(in(gb), fs)
            end
        end
    end
end
