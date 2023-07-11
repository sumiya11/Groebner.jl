
const loglevel = 0

@testset "Learn & apply" begin
    K = AbstractAlgebra.GF(2^31-1)
    R, (x, y) = PolynomialRing(K, ["x", "y"], ordering=:degrevlex)
    R2, xs = PolynomialRing(K, ["x$i" for i in 1:30], ordering=:degrevlex)

    @test_throws DomainError Groebner.groebner_learn([R(0), R(0)], loglevel=loglevel)
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
        (system=[sum(xs) + prod(xs), sum(xs)^2, prod(xs) - 1],),
    ]
    for case in cases
        # Learn and apply on the same system
        system = case.system
        true_gb = Groebner.groebner(system, loglevel=loglevel)
        graph, gb_1 = Groebner.groebner_learn(system, loglevel=loglevel)
        flag, gb_2 = Groebner.groebner_apply!(graph, system, loglevel=loglevel)
        @test flag && gb_2 == true_gb

        # Apply on a different system N times!
        N = 5
        X = gens(parent(first(system)))
        for _ in 1:N
            point = map(t -> iszero(t) ? t + 1 : t, rand(K, length(X))) .* X
            system_ = map(f -> evaluate(f, point), system)        
            true_gb = Groebner.groebner(system_, loglevel=loglevel)
            flag, gb_2 = Groebner.groebner_apply!(graph, system_, loglevel=loglevel)
            @test flag && gb_2 == true_gb
        end
    end
end

@testset "Learn & apply tricky" begin
    K = AbstractAlgebra.GF(2^31-1)
    R, (x, y) = PolynomialRing(K, ["x", "y"], ordering=:degrevlex)
    
    # s-poly of x + 1 and x*y + 7y is y - 7y.
    system_1 = [x + 1, x*y + 7y]
    graph, gb_1 = Groebner.groebner_learn(system_1, loglevel=loglevel)
    flag, gb_2 = Groebner.groebner_apply!(graph, system_1, loglevel=loglevel)
    @test flag && gb_2 == gb_1

    # TODO: Think about something here
    # @test_throws AssertionError Groebner.groebner_apply!(graph, [x, y], loglevel=loglevel)

    # Cancellation of a leading term:
    # s-poly of x + 1 and x*y + y is 0 = y - y.
    system_2 = [x + 1, x*y + y]
    flag, gb_2 = Groebner.groebner_apply!(graph, system_2, loglevel=loglevel)
    @test !flag

    # s-poly of x + y + 1 and x*y + 7y is y^2 - 7y
    system_1 = [x + y + 1, x*y + 7y]
    graph, gb_1 = Groebner.groebner_learn(system_1, loglevel=loglevel)
    flag, gb_2 = Groebner.groebner_apply!(graph, system_1, loglevel=loglevel)
    @test flag && gb_2 == gb_1

    # Cancellation of a trailing term:
    # s-poly of x + y + 1 and x*y + y is y^2
    # TODO: produce a warning here
    system_2 = [x + y + 1, x*y + y]
    flag, gb_2 = Groebner.groebner_apply!(graph, system_2, loglevel=loglevel)
    @test !flag

    # Input is a Groebner basis already:
    N = 3
    for system in [
        Groebner.noonn(4, ordering=:degrevlex, ground=K), 
        Groebner.katsuran(5, ordering=:degrevlex, ground=K),
        [x, x^2, y^2, x*y, x^3, y^4, x^10, y^10, x*y^10]
    ]
        gb = Groebner.groebner(system, loglevel=loglevel)
        for i in 1:N
            graph, gb_1 = Groebner.groebner_learn(gb, loglevel=loglevel)
            flag, gb_2 = Groebner.groebner_apply!(graph, gb, loglevel=loglevel)
            @test flag && gb_2 == gb_1 == gb
        end
    end 
end
