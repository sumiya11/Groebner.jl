import Random

params = (loglevel=0, sweep=true)
# TODO: do not turn this off 
Groebner.invariants_enabled() = true

@testset "Learn & apply" begin
    K = AbstractAlgebra.GF(2^31 - 1)
    K2 = AbstractAlgebra.GF(2^30 + 3)

    R, (x, y) = PolynomialRing(K, ["x", "y"], ordering=:degrevlex)
    R2, (x2, y2) = PolynomialRing(K2, ["x", "y"], ordering=:degrevlex)
    R2, xs = PolynomialRing(K, ["x$i" for i in 1:30], ordering=:degrevlex)

    @test_throws DomainError Groebner.groebner_learn([R(0), R(0)]; params...)
    cases = [
        (system=[x, R(0)],),
        (system=[R(1)],),
        (system=[x, y],),
        (system=[x, x],),
        (system=[x, y, x * y],),
        (system=[x + 1, y + 2, x * y + 3],),
        (system=[x^20 * y + x + 1, x * y^20 + y + 1],),
        (system=[x^20 * y + x + 1, x * y^20 + y + 1],),
        (system=Groebner.noonn(3, ordering=:degrevlex, ground=K),),
        (system=Groebner.noonn(4, ordering=:degrevlex, ground=K),),
        (system=Groebner.noonn(5, ordering=:degrevlex, ground=K),),
        (system=Groebner.katsuran(3, ordering=:degrevlex, ground=K),),
        (system=Groebner.katsuran(4, ordering=:degrevlex, ground=K),),
        (system=Groebner.katsuran(5, ordering=:degrevlex, ground=K),),
        (system=Groebner.cyclicn(5, ordering=:degrevlex, ground=K),),
        (system=Groebner.rootn(5, ordering=:lex, ground=K),),
        (system=Groebner.rootn(5, ordering=:deglex, ground=K),),
        (system=Groebner.rootn(6, ordering=:degrevlex, ground=K),),
        (system=Groebner.eco5(ordering=:degrevlex, ground=K),),
        (system=Groebner.ku10(ordering=:degrevlex, ground=K),),
        (system=Groebner.kinema(ordering=:degrevlex, ground=K),),
        (system=Groebner.sparse5(ordering=:degrevlex, ground=K),),
        (system=Groebner.s9_1(ordering=:degrevlex, ground=K),),
        (system=[sum(xs) + prod(xs), sum(xs)^2, prod(xs) - 1],)
    ]
    for case in cases
        # Learn and apply on the same system
        system = case.system
        true_gb = Groebner.groebner(system; params...)
        graph, gb_1 = Groebner.groebner_learn(system; params...)
        flag, gb_2 = Groebner.groebner_apply!(graph, system; params...)
        @test flag && gb_2 == true_gb

        # Apply on a different system N times
        N = 5
        X = gens(parent(first(system)))
        for _ in 1:N
            point = map(t -> iszero(t) ? t + 1 : t, rand(K, length(X))) .* X
            system_ = map(f -> evaluate(f, point), system)
            true_gb = Groebner.groebner(system_; params...)
            flag, gb_2 = Groebner.groebner_apply!(graph, system_; params...)
            @test flag && gb_2 == true_gb
        end

        # Apply on the same system but modulo a different prime 
        system_2 = map(f -> map_coefficients(c -> K2(data(c)), f), system)
        try
            flag, gb_2 = Groebner.groebner_apply!(graph, system_2; params...)
            @test flag && gb_2 == Groebner.groebner(system_2)
        catch error
            if error isa AssertionError
                @test_broken false
            else
                rethrow(error)
            end
        end
    end
end

@testset "Learn & apply, orderings" begin
    K = GF(2^31 - 1)
    R, (x, y) = PolynomialRing(K, ["x", "y"], ordering=:lex)

    ord_1 = Groebner.Lex()
    graph_1, gb_1 = Groebner.groebner_learn([x + 2y + 3, y], ordering=ord_1)
    @test gb_1 == [y, x + 3]

    flag, gb_1_apply = Groebner.groebner_apply!(graph_1, [x + y + 3, y])
    @test flag && gb_1 == gb_1_apply

    ord_2 = Groebner.Lex(y, x)
    graph_2, gb_2 = Groebner.groebner_learn([y, x + 2y^2 + 3], ordering=ord_2)
    @test gb_2 == [x + 3, y]

    flag, gb_2_apply = Groebner.groebner_apply!(graph_2, [y, x + 2y^2 + 3])
    @test flag && gb_2 == gb_2_apply

    K = GF(2^31 - 1)
    n = 10
    R, x = PolynomialRing(K, [["x$i" for i in 1:n]...], ordering=:degrevlex)
    F = (x .+ (1:n) .* circshift(x, 1)) .^ 2
    # F = [(x1 + xn)^2, (x2 + 2 x1)^2, ..., (xn + n x_{n-1})^2]
    for i in 0:(2n)
        if i < n
            xi = circshift(x, -i)
        else
            xi = Random.shuffle(x)
        end

        ord = Groebner.Lex(xi)
        input = Random.shuffle(F)
        gb = Groebner.groebner(input, ordering=ord)
        graph, gb_1 = Groebner.groebner_learn(input, ordering=ord)
        flag, gb_2 = Groebner.groebner_apply!(graph, input)
        @test gb == gb_1
        @test flag && gb == gb_2

        ord = Groebner.DegRevLex(xi)
        input = Random.shuffle(F)
        gb = Groebner.groebner(input, ordering=ord)
        graph, gb_1 = Groebner.groebner_learn(input, ordering=ord)
        flag, gb_2 = Groebner.groebner_apply!(graph, input)
        @test gb == gb_1
        @test flag && gb == gb_2
    end
end

@testset "Learn & apply, tricky" begin
    for K in [GF(2^31 - 1), GF(2^62 + 135)]
        R, (x, y) = PolynomialRing(K, ["x", "y"], ordering=:degrevlex)

        s = [x^100 * y + y^100, x * y^100 + y]
        graph, gb_1 = Groebner.groebner_learn(s; params...)
        flag, gb_2 = Groebner.groebner_apply!(graph, s; params...)
        @test gb_1 == gb_2 == [x * y^100 + y, x^100 * y + y^100, y^199 - x^99 * y]

        s = [x^200 * y + y^200, x * y^200 + y]
        graph, gb_1 = Groebner.groebner_learn(s; params...)
        flag, gb_2 = Groebner.groebner_apply!(graph, s; params...)
        @test gb_1 == gb_2 == [x * y^200 + y, x^200 * y + y^200, y^399 - x^199 * y]

        s = [x^1000 * y + y^1000, x * y^1000 + y]
        graph, gb_1 = Groebner.groebner_learn(s; params...)
        for _ in 1:5
            flag, gb_2 = Groebner.groebner_apply!(graph, s; params...)
            @test gb_1 == gb_2 == [x * y^1000 + y, x^1000 * y + y^1000, y^1999 - x^999 * y]
        end

        # TODO
        # @test_throws AssertionError Groebner.groebner_apply!(graph, s; sweep=!params.sweep)
    end

    K = AbstractAlgebra.GF(2^31 - 1)
    R, (x, y) = PolynomialRing(K, ["x", "y"], ordering=:degrevlex)

    # s-poly of x + 1 and x*y + 7y is y - 7y.
    system_1 = [x + 1, x * y + 7y]
    graph, gb_1 = Groebner.groebner_learn(system_1; params...)
    flag, gb_2 = Groebner.groebner_apply!(graph, system_1; params...)
    @test flag && gb_2 == gb_1

    # TODO: Think about something here
    # @test_throws AssertionError Groebner.groebner_apply!(graph, [x, y]; params...)

    # Cancellation of a leading term:
    # s-poly of x + 1 and x*y + y is 0 = y - y.
    system_2 = [x + 1, x * y + y]
    flag, gb_2 = Groebner.groebner_apply!(graph, system_2; params...)
    @test !flag

    # s-poly of x + y + 1 and x*y + 7y is y^2 - 7y
    system_1 = [x + y + 1, x * y + 7y]
    graph, gb_1 = Groebner.groebner_learn(system_1; params...)
    flag, gb_2 = Groebner.groebner_apply!(graph, system_1; params...)
    @test flag && gb_2 == gb_1

    # Cancellation of a trailing term:
    # s-poly of x + y + 1 and x*y + y is y^2
    # TODO: produce a warning here
    system_2 = [x + y + 1, x * y + y]
    flag, gb_2 = Groebner.groebner_apply!(graph, system_2; params...)
    @test !flag

    # Input is a Groebner basis already:
    N = 3
    for system in [
        Groebner.noonn(4, ordering=:degrevlex, ground=K),
        Groebner.katsuran(5, ordering=:degrevlex, ground=K),
        [x, x^2, y^2, x * y, x^3, y^4, x^10, y^10, x * y^10]
    ]
        gb = Groebner.groebner(system; params...)
        for i in 1:N
            graph, gb_1 = Groebner.groebner_learn(gb; params...)
            flag, gb_2 = Groebner.groebner_apply!(graph, gb; params...)
            @test flag && gb_2 == gb_1 == gb
        end
    end
end
