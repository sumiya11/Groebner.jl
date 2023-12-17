import Random

params = (loglevel=0, sweep=true)

@testset "learn & apply, same field" begin
    K = AbstractAlgebra.GF(2^31 - 1)
    R, (x, y) = polynomial_ring(K, ["x", "y"])
    trace, gb1 = Groebner.groebner_learn([x, y])
    flag, gb2 = Groebner.groebner_apply!(trace, [x, y])
    @assert gb1 == gb2 && flag

    K = AbstractAlgebra.GF(2^31 - 1)
    R, (x, y) = polynomial_ring(K, ["x", "y"], ordering=:degrevlex)
    R2, xs = polynomial_ring(K, ["x$i" for i in 1:30], ordering=:degrevlex)

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
        trace, gb_1 = Groebner.groebner_learn(system; params...)
        flag, gb_2 = Groebner.groebner_apply!(trace, system; params...)
        @test flag && gb_2 == true_gb

        # Apply on a different system N times
        N = 5
        X = gens(parent(first(system)))
        for _ in 1:N
            point = map(t -> iszero(t) ? t + 1 : t, rand(K, length(X))) .* X
            system_ = map(f -> evaluate(f, point), system)
            true_gb = Groebner.groebner(system_; params...)
            flag, gb_2 = Groebner.groebner_apply!(trace, system_; params...)
            @test flag && gb_2 == true_gb
        end
    end
end

@testset "learn & apply, different modulo" begin
    # Some small tests and corner cases
    K, K2, K3 = GF(5), GF(7), GF(11)
    R, (x, y) = polynomial_ring(K, ["x", "y"], ordering=:degrevlex)
    R2, (x2, y2) = polynomial_ring(K2, ["x", "y"], ordering=:degrevlex)
    R3, (x3, y3) = polynomial_ring(K3, ["x", "y"], ordering=:degrevlex)
    system = [x^4 + 2y^3 + 3x^2 + 4y, 4y^4 + 3x^3 + 2y^2 + 1x]
    system2 = map(f -> map_coefficients(c -> K2(data(c)), f), system)

    trace, gb_1 = Groebner.groebner_learn(system; params...)
    flag, gb_2 = Groebner.groebner_apply!(trace, system2; params...)
    @test flag && gb_2 == Groebner.groebner(system2)

    trace, gb_1 = Groebner.groebner_learn(system2; params...)
    flag, gb_2 = Groebner.groebner_apply!(trace, system; params...)
    @test flag && gb_2 == Groebner.groebner(system)

    system = [6x2^4 + 6y2 + 1, y2 + 3]
    system2 = map(f -> map_coefficients(c -> K(data(c)), f), system)
    system3 = map(f -> map_coefficients(c -> K3(data(c)), f), system)
    trace, gb_1 = Groebner.groebner_learn(system; params...)
    flag, gb_2 = Groebner.groebner_apply!(trace, system2; params...)
    @test flag && gb_2 == Groebner.groebner(system2)
    flag, gb_3 = Groebner.groebner_apply!(trace, system3; params...)
    @test flag && gb_3 == Groebner.groebner(system3)

    system = [5x2 + 5, 5y2 + 1]
    system2 = map(f -> map_coefficients(c -> K(data(c)), f), system)
    trace, gb_1 = Groebner.groebner_learn(system; params...)
    @test_throws DomainError Groebner.groebner_apply!(trace, system2; params...)
    system = [x2 + 5, y2 + 1]
    system2 = map(f -> map_coefficients(c -> K(data(c)), f), system)
    trace, gb_1 = Groebner.groebner_learn(system; params...)
    @test_throws DomainError Groebner.groebner_apply!(trace, system2; params...)

    # Going from small characteristic to large is not allowed
    # NOTE: it should be allowed
    K1, K2 = GF(2^31 - 1), GF(2^60 + 33)
    R, (x, y) = polynomial_ring(K1, ["x", "y"], ordering=:degrevlex)
    R2, (x2, y2) = polynomial_ring(K2, ["x", "y"], ordering=:degrevlex)
    system2 = [(2^49 + 1) * x2 - 1, (2^50) * y2 + (2^56 + 99)]
    system = map(f -> map_coefficients(c -> K1(data(c)), f), system)
    trace, gb_1 = Groebner.groebner_learn(system; params...)
    flag, gb_2 = Groebner.groebner_apply!(trace, system2; params...)
    @test gb_2 == [y2 + (2^56 + 99) // K2(2^50), x2 - 1 // K2(2^49 + 1)]

    # Going from one monomial ordering to another is not allowed
    K1, K2 = GF(2^31 - 1), GF(2^60 + 33)
    R, (x, y) = polynomial_ring(K1, ["x", "y"], ordering=:lex)
    R2, (x2, y2) = polynomial_ring(K2, ["x", "y"], ordering=:degrevlex)
    system = [x + 1, y - 1]
    system2 = [x2 + 1, y2 - 1]
    trace, gb_1 = Groebner.groebner_learn(system; params...)
    @test_throws DomainError Groebner.groebner_apply!(trace, system2; params...)

    K1, K2 = GF(2^31 - 1), GF(2^60 + 33)
    R, (x, y) = polynomial_ring(K1, ["x", "y"], ordering=:lex)
    R2, (x2, y2) = polynomial_ring(K2, ["x", "y"], ordering=:degrevlex)
    system = [x + 1, y - 1]
    system2 = [x2 + 1, y2 - 1]
    trace, gb_1 = Groebner.groebner_learn(system; params...)
    @test_throws AssertionError Groebner.groebner_apply!(
        trace,
        system2;
        ordering=Groebner.DegRevLex(),
        params...
    )

    # NOTE: these systems are only relatively very large.
    # We should also test for larger systems and larger primes!!
    R, (x, y) = polynomial_ring(ZZ, ["x", "y"], ordering=:degrevlex)
    cases = [
        (system=[x, R(0)],),
        (system=[R(1)],),
        (system=[x, y],),
        (system=[x, x],),
        (system=[x, y, x * y],),
        (system=[x + 1, y + 2, x * y + 3],),
        (system=[x^20 * y + x + 1, x * y^20 + y + 1],),
        (system=[x^20 * y + x + 1, x * y^20 + y + 1],),
        (system=Groebner.noonn(3, ordering=:lex, ground=ZZ),),
        (system=Groebner.noonn(4, ordering=:deglex, ground=ZZ),),
        (system=Groebner.noonn(5, ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.katsuran(3, ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.katsuran(4, ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.katsuran(5, ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.katsuran(6, ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.katsuran(7, ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.cyclicn(5, ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.cyclicn(6, ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.rootn(5, ordering=:lex, ground=ZZ),),
        (system=Groebner.rootn(5, ordering=:deglex, ground=ZZ),),
        (system=Groebner.rootn(6, ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.eco5(ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.ku10(ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.kinema(ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.sparse5(ordering=:degrevlex, ground=ZZ),),
        (system=Groebner.s9_1(ordering=:degrevlex, ground=ZZ),)
    ]

    # Some bigger tests
    K = AbstractAlgebra.GF(2^31 - 1)
    Ks = [
        AbstractAlgebra.GF(2^31 - 1),
        AbstractAlgebra.GF(2^30 + 3),
        AbstractAlgebra.GF(2^31 + 11),
        AbstractAlgebra.GF(2^20 + 7),
        AbstractAlgebra.GF(2^20 + 13),
        AbstractAlgebra.GF(2^20 + 25),
        AbstractAlgebra.GF(2^15 + 3),
        AbstractAlgebra.GF(2^15 + 11),
        AbstractAlgebra.GF(2^15 + 15)
        # it may be dangerous to use learn & apply for small primes..
        # AbstractAlgebra.GF(2^11 + 5),
        # AbstractAlgebra.GF(2^11 + 15),
        # AbstractAlgebra.GF(2^11 + 21)
    ]
    for case in cases
        # Learn and apply on the same system
        system = case.system
        system_zp = map(f -> map_coefficients(c -> K(BigInt(c)), f), system)
        true_gb = Groebner.groebner(system_zp; params...)
        trace, gb_1 = Groebner.groebner_learn(system_zp; params...)
        flag, gb_2 = Groebner.groebner_apply!(trace, system_zp; params...)
        @test flag && gb_2 == true_gb

        # Apply on the same system but modulo a different prime 
        for K_ in Ks
            system_zp_2 = map(f -> map_coefficients(c -> K_(BigInt(c)), f), system)
            flag, gb_2 = Groebner.groebner_apply!(trace, system_zp_2; params...)
            @test flag && gb_2 == Groebner.groebner(system_zp_2)
        end
    end

    K1 = GF(2^31 - 1)
    R, (x, y) = polynomial_ring(ZZ, ["x", "y"], ordering=:lex)
    system = [BigInt(2)^1000 * x + (BigInt(2)^1001 + 1) * y + 1]
    system_zp = map(f -> map_coefficients(c -> K1(BigInt(c)), f), system)
    trace, gb_1 =
        Groebner.groebner_learn(system_zp; ordering=Groebner.DegRevLex(y, x), params...)
    for K in Primes.nextprimes(2^31, 200)
        system_zp = map(f -> map_coefficients(c -> GF(K)(BigInt(c)), f), system)
        true_gb = Groebner.groebner(system_zp; ordering=Groebner.DegRevLex(y, x), params...)
        flag, gb_2 = Groebner.groebner_apply!(trace, system_zp; params...)
        @test flag && gb_2 == true_gb
    end
end

@testset "learn & apply, orderings" begin
    K = GF(2^31 - 1)
    R, (x, y) = polynomial_ring(K, ["x", "y"], ordering=:lex)

    ord_1 = Groebner.Lex()
    trace_1, gb_1 = Groebner.groebner_learn([x + 2y + 3, y], ordering=ord_1)
    @test gb_1 == [y, x + 3]

    flag, gb_1_apply = Groebner.groebner_apply!(trace_1, [x + y + 3, y])
    @test flag && gb_1 == gb_1_apply

    ord_2 = Groebner.Lex(y, x)
    trace_2, gb_2 = Groebner.groebner_learn([y, x + 2y^2 + 3], ordering=ord_2)
    @test gb_2 == [x + 3, y]

    flag, gb_2_apply = Groebner.groebner_apply!(trace_2, [y, x + 2y^2 + 3])
    @test flag && gb_2 == gb_2_apply

    K = GF(2^31 - 1)
    n = 10
    R, x = polynomial_ring(K, [["x$i" for i in 1:n]...], ordering=:degrevlex)
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
        trace, gb_1 = Groebner.groebner_learn(input, ordering=ord)
        flag, gb_2 = Groebner.groebner_apply!(trace, input)
        @test gb == gb_1
        @test flag && gb == gb_2

        ord = Groebner.DegRevLex(xi)
        input = Random.shuffle(F)
        gb = Groebner.groebner(input, ordering=ord)
        trace, gb_1 = Groebner.groebner_learn(input, ordering=ord)
        flag, gb_2 = Groebner.groebner_apply!(trace, input)
        @test gb == gb_1
        @test flag && gb == gb_2
    end
end

@testset "learn & apply, tricky" begin
    for K in [GF(2^31 - 1), GF(2^62 + 135)]
        R, (x, y) = polynomial_ring(K, ["x", "y"], ordering=:degrevlex)

        s = [x^100 * y + y^100, x * y^100 + y]
        trace, gb_1 = Groebner.groebner_learn(s; params...)
        flag, gb_2 = Groebner.groebner_apply!(trace, s; params...)
        @test gb_1 == gb_2 == [x * y^100 + y, x^100 * y + y^100, y^199 - x^99 * y]

        s = [x^200 * y + y^200, x * y^200 + y]
        trace, gb_1 = Groebner.groebner_learn(s; params...)
        flag, gb_2 = Groebner.groebner_apply!(trace, s; params...)
        @test gb_1 == gb_2 == [x * y^200 + y, x^200 * y + y^200, y^399 - x^199 * y]

        s = [x^1000 * y + y^1000, x * y^1000 + y]
        trace, gb_1 = Groebner.groebner_learn(s; params...)
        for _ in 1:5
            flag, gb_2 = Groebner.groebner_apply!(trace, s; params...)
            @test gb_1 == gb_2 == [x * y^1000 + y, x^1000 * y + y^1000, y^1999 - x^999 * y]
        end

        @test_throws DomainError Groebner.groebner_apply!(trace, s; sweep=!params.sweep)
    end

    K = AbstractAlgebra.GF(2^31 - 1)
    R, (x, y) = polynomial_ring(K, ["x", "y"], ordering=:degrevlex)

    # s-poly of x + 1 and x*y + 7y is y - 7y.
    system_1 = [x + 1, x * y + 7y]
    trace, gb_1 = Groebner.groebner_learn(system_1; params...)
    flag, gb_2 = Groebner.groebner_apply!(trace, system_1; params...)
    @test flag && gb_2 == gb_1

    @test_throws DomainError Groebner.groebner_apply!(trace, [x, y]; params...)

    # Cancellation of a leading term:
    # s-poly of x + 1 and x*y + y is 0 = y - y.
    system_2 = [x + 1, x * y + y]
    flag, gb_2 = Groebner.groebner_apply!(trace, system_2; params...)
    @test !flag

    # s-poly of x + y + 1 and x*y + 7y is y^2 - 7y
    system_1 = [x + y + 1, x * y + 7y]
    trace, gb_1 = Groebner.groebner_learn(system_1; params...)
    flag, gb_2 = Groebner.groebner_apply!(trace, system_1; params...)
    @test flag && gb_2 == gb_1

    # Cancellation of a trailing term:
    # s-poly of x + y + 1 and x*y + y is y^2
    # TODO: produce a warning here
    system_2 = [x + y + 1, x * y + y]
    flag, gb_2 = Groebner.groebner_apply!(trace, system_2; params...)
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
            trace, gb_1 = Groebner.groebner_learn(gb; params...)
            flag, gb_2 = Groebner.groebner_apply!(trace, gb; params...)
            @test flag && gb_2 == gb_1 == gb
        end
    end
end
