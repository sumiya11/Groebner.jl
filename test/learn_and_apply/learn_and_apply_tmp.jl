import Random

params = (loglevel=0, sweep=true)

@testset "learn & apply, same field" begin
    K = AbstractAlgebra.GF(2^31 - 1)
    R, (x, y) = polynomial_ring(K, ["x", "y"], ordering=:degrevlex)
    trace, gb1 = Groebner.groebner_learn([x + 2, y + 3])
    flag, cfs = Groebner.groebner_applyX!(
        trace,
        [[UInt32(1), UInt32(3)], [UInt32(1), UInt32(2)]],
        UInt32(2^27 + 39)
    )
    @assert flag && cfs == [[UInt32(1), UInt32(2)], [UInt32(1), UInt32(3)]]

    cases = [
        (system=[x],),
        (system=[R(1)],),
        (system=[x, y],),
        (system=[x, x],),
        (system=[x, y, x * y],),
        (system=[x + 1, y + 2, x * y + 3],),
        (system=[x^20 * y + x + 1, x * y^20 + y + 1],),
        (system=[x^20 * y + x + 1, x * y^20 + y + 1],),
        (system=Groebner.noonn(3, ordering=:degrevlex, k=K),),
        (system=Groebner.noonn(4, ordering=:degrevlex, k=K),),
        (system=Groebner.noonn(5, ordering=:degrevlex, k=K),),
        (system=Groebner.katsuran(3, ordering=:degrevlex, k=K),),
        (system=Groebner.katsuran(4, ordering=:degrevlex, k=K),),
        (system=Groebner.katsuran(5, ordering=:degrevlex, k=K),),
        (system=Groebner.cyclicn(5, ordering=:degrevlex, k=K),),
        (system=Groebner.rootn(5, ordering=:degrevlex, k=K),),
        (system=Groebner.rootn(5, ordering=:degrevlex, k=K),),
        (system=Groebner.rootn(6, ordering=:degrevlex, k=K),),
        (system=Groebner.eco5(ordering=:degrevlex, k=K),),
        (system=Groebner.eco10(ordering=:degrevlex, k=K),),
        (system=Groebner.ku10(ordering=:degrevlex, k=K),),
        (system=Groebner.kinema(ordering=:degrevlex, k=K),),
        (system=Groebner.sparse5(ordering=:degrevlex, k=K),),
        (system=Groebner.s9_1(ordering=:degrevlex, k=K),)
    ]

    for case in cases
        # Learn and apply on the same system
        system = case.system
        true_gb = Groebner.groebner(system; params...)
        trace, gb_1 = Groebner.groebner_learn(system; params...)
        cfs = map(
            f -> map(c -> UInt32(data(c)), collect(AbstractAlgebra.coefficients(f))),
            system
        )
        flag, cfs_2 = Groebner.groebner_applyX!(
            trace,
            cfs,
            UInt32(AbstractAlgebra.characteristic(K));
            params...
        )
        @test flag && map(f -> collect(coefficients(f)), true_gb) == cfs_2
    end
end

@testset "learn & apply, different field" begin
    # Some small tests and corner cases
    K, K2, K3 = GF(5), GF(7), GF(11)
    R, (x, y) = polynomial_ring(K, ["x", "y"], ordering=:degrevlex)
    R2, (x2, y2) = polynomial_ring(K2, ["x", "y"], ordering=:degrevlex)
    R3, (x3, y3) = polynomial_ring(K3, ["x", "y"], ordering=:degrevlex)
    system = [x^4 + 2y^3 + 3x^2 + 4y, 4y^4 + 3x^3 + 2y^2 + 1x]
    system2 = map(f -> map_coefficients(c -> K2(data(c)), f), system)

    trace, gb_1 = Groebner.groebner_learn(system; params...)
    cfs = map(
        f -> map(c -> UInt32(data(c)), collect(AbstractAlgebra.coefficients(f))),
        system2
    )
    flag, cfs_2 = Groebner.groebner_applyX!(
        trace,
        cfs,
        UInt32(AbstractAlgebra.characteristic(K2));
        params...
    )
    @test flag && map(f -> collect(coefficients(f)), Groebner.groebner(system2)) == cfs_2

    # NOTE: these systems are only relatively very large.
    # We should also test for larger systems and larger primes!!
    R, (x, y) = polynomial_ring(ZZ, ["x", "y"], ordering=:degrevlex)
    cases = [
        (system=[x],),
        (system=[R(1)],),
        (system=[x, y],),
        (system=[x, x],),
        (system=[x, y, x * y],),
        (system=[x + 1, y + 2, x * y + 3],),
        (system=[x^20 * y + x + 1, x * y^20 + y + 1],),
        (system=[x^20 * y + x + 1, x * y^20 + y + 1],),
        (system=Groebner.noonn(3, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.noonn(4, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.noonn(5, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.katsuran(3, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.katsuran(4, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.katsuran(5, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.katsuran(6, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.katsuran(7, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.cyclicn(5, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.cyclicn(6, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.rootn(5, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.rootn(5, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.rootn(6, ordering=:degrevlex, k=ZZ),),
        (system=Groebner.eco5(ordering=:degrevlex, k=ZZ),),
        (system=Groebner.eco10(ordering=:degrevlex, k=ZZ),),
        (system=Groebner.eco11(ordering=:degrevlex, k=ZZ),),
        (system=Groebner.ku10(ordering=:degrevlex, k=ZZ),),
        (system=Groebner.kinema(ordering=:degrevlex, k=ZZ),),
        (system=Groebner.sparse5(ordering=:degrevlex, k=ZZ),),
        (system=Groebner.s9_1(ordering=:degrevlex, k=ZZ),)
    ]

    # Some bigger tests
    K = AbstractAlgebra.GF(2^31 - 1)
    Ks = [
        AbstractAlgebra.GF(2^31 - 1),
        AbstractAlgebra.GF(2^30 + 3),
        AbstractAlgebra.GF(2^31 + 11),
        AbstractAlgebra.GF(2^27 - 39),
        AbstractAlgebra.GF(2^20 + 7),
        AbstractAlgebra.GF(2^20 + 13),
        AbstractAlgebra.GF(2^20 + 25),
    ]
    for case in cases
        # Learn and apply on the same system
        system = case.system
        system_zp = map(f -> map_coefficients(c -> K(BigInt(c)), f), system)
        true_gb = Groebner.groebner(system_zp; params...)
        trace, gb_1 = Groebner.groebner_learn(system_zp; params...)

        # Apply on the same system but modulo a different prime 
        for K_ in Ks
            system_zp_2 = map(f -> map_coefficients(c -> K_(BigInt(c)), f), system)
            cfs = map(
                f ->
                    map(c -> UInt32(data(c)), collect(AbstractAlgebra.coefficients(f))),
                system_zp_2
            )
            flag, cfs_2 = Groebner.groebner_applyX!(
                trace,
                cfs,
                UInt32(AbstractAlgebra.characteristic(K_));
                params...
            )
            @test flag &&
                  map(f -> collect(coefficients(f)), Groebner.groebner(system_zp_2)) ==
                  cfs_2
        end
    end
end
