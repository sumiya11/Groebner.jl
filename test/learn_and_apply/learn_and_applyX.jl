import Random

params = (loglevel=0, sweep=true)

@testset "learn & apply, same field" begin
    K = AbstractAlgebra.GF(2^31 - 1)
    R, (x, y) = polynomial_ring(K, ["x", "y"], internal_ordering=:degrevlex)
    trace, gb1 = Groebner.groebner_learn([x + 2, y + 3])
    flag, cfs = Groebner.groebner_applyX!(
        trace,
        [[UInt32(1), UInt32(3)], [UInt32(1), UInt32(2)]],
        UInt32(2^27 - 39)
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
        (system=Groebner.Examples.noonn(3, internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.noonn(4, internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.noonn(5, internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.katsuran(3, internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.katsuran(4, internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.katsuran(5, internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.cyclicn(5, internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.rootn(5, internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.rootn(5, internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.rootn(6, internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.eco5(internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.eco10(internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.ku10(internal_ordering=:degrevlex, k=K),),
        (system=Groebner.Examples.kinema(internal_ordering=:degrevlex, k=K),),
        (system=Groebner.s9_1(internal_ordering=:degrevlex, k=K),)
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
    K, K2, K3 = GF(2^30 + 3), GF(2^31 - 1), GF(2^27 - 39)
    R, (x, y) = polynomial_ring(K, ["x", "y"], internal_ordering=:degrevlex)
    R2, (x2, y2) = polynomial_ring(K2, ["x", "y"], internal_ordering=:degrevlex)
    R3, (x3, y3) = polynomial_ring(K3, ["x", "y"], internal_ordering=:degrevlex)
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
    R, (x, y) = polynomial_ring(ZZ, ["x", "y"], internal_ordering=:degrevlex)
    cases = [
        (system=[x],),
        (system=[R(1)],),
        (system=[x, y],),
        (system=[x, x],),
        (system=[x, y, x * y],),
        (system=[x + 1, y + 2, x * y + 3],),
        (system=[x^20 * y + x + 1, x * y^20 + y + 1],),
        (system=[x^20 * y + x + 1, x * y^20 + y + 1],),
        (system=Groebner.Examples.noonn(3, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.noonn(4, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.noonn(5, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.katsuran(3, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.katsuran(4, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.katsuran(5, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.katsuran(6, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.katsuran(7, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.cyclicn(5, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.cyclicn(6, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.rootn(5, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.rootn(5, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.rootn(6, internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.eco5(internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.eco10(internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.eco11(internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.ku10(internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.Examples.kinema(internal_ordering=:degrevlex, k=ZZ),),
        (system=Groebner.s9_1(internal_ordering=:degrevlex, k=ZZ),)
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
        AbstractAlgebra.GF(2^20 + 25)
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

    # Leading term in the matrix cancel out
    R, (x, y) = polynomial_ring(QQ, ["x", "y"], internal_ordering=:degrevlex)
    p, p2 = 3, 19
    sys = [-11 * x * y + 53 * y, 83 * x * y + x - 70 * y]
    sys_mod_p = map(
        f -> map_coefficients(c -> GF(p)(numerator(c)), f),
        [x * y + 2 * y, 2 * x * y + x + 2 * y]
    )
    sys_mod_p2 = map(
        f -> map_coefficients(c -> GF(p2)(numerator(c)), f),
        [8 * x * y + 15 * y, 7 * x * y + x + 6 * y]
    )
    Groebner.groebner(sys_mod_p)
    Groebner.groebner(sys_mod_p2)

    _cfs(s) = map(f -> map(UInt32 ∘ data, collect(AbstractAlgebra.coefficients(f))), s)

    trace, gb1 = Groebner.groebner_learn(sys_mod_p)
    flag2, gb2 = Groebner.groebner_applyX!(trace, _cfs(sys_mod_p), UInt32(3))
    @test flag
    flag3, gb3 = Groebner.groebner_applyX!(trace, _cfs(sys_mod_p2), UInt32(19))
    flag4, gb4 = Groebner.groebner_applyX!(trace, _cfs(sys_mod_p2), UInt32(19))
    flag5, gb5 = Groebner.groebner_applyX!(trace, _cfs(sys_mod_p2), UInt32(19))

    trace, gb = Groebner.groebner_learn(sys_mod_p2)
    flag, gb = Groebner.groebner_applyX!(trace, _cfs(sys_mod_p2), UInt32(19))
    @test flag
    flag, gb = Groebner.groebner_applyX!(trace, _cfs(sys_mod_p), UInt32(3))
end
