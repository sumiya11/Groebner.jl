import AbstractAlgebra, Groebner

@info "Loading benchmark models.."

include((@__DIR__) * "/generate/benchmark_systems/SIAN/SIAN.jl")
include((@__DIR__) * "/generate/benchmark_systems/MQ/MQ.jl")
include((@__DIR__) * "/generate/benchmark_systems/SI/SI.jl")
include((@__DIR__) * "/generate/benchmark_systems/biomodels/BIOMD.jl")

function get_benchmark_suite(id)
    s = Symbol(:benchmark_set_, Symbol(id))
    @eval $s()
end

function dummy_system(name, ground_field)
    ring, (x, y) = AbstractAlgebra.polynomial_ring(ground_field, ["x", "y"])
    (name, [x^2 + y^2 - 1, x + y])
end

function benchmark_set_0()
    (
        name="dummy benchmark set",
        field=AbstractAlgebra.GF(7),
        systems=[
            dummy_system("dummy 1", AbstractAlgebra.GF(7)),
            dummy_system("dummy 2", AbstractAlgebra.GF(7))
        ]
    )
end

function benchmark_set_1()
    ground_field = AbstractAlgebra.GF(2^30 + 3)

    systems = [
        dummy_system("dummy", ground_field),
        ("cyclic 7", Groebner.cyclicn(7, k=ground_field)),
        ("cyclic 8", Groebner.cyclicn(8, k=ground_field)),
        ("cyclic 9", Groebner.cyclicn(9, k=ground_field)),
        # ("cyclic 10", Groebner.cyclicn(10, k=ground_field)),
        ("katsura 10", Groebner.katsuran(10, k=ground_field)),
        ("katsura 11", Groebner.katsuran(11, k=ground_field)),
        ("katsura 12", Groebner.katsuran(12, k=ground_field)),
        ("katsura 13", Groebner.katsuran(13, k=ground_field)),
        # ("katsura 14", Groebner.katsuran(14, k=ground_field)),
        ("eco 11", Groebner.eco11(k=ground_field)),
        ("eco 12", Groebner.eco12(k=ground_field)),
        ("eco 13", Groebner.eco13(k=ground_field)),
        ("eco 14", Groebner.eco14(k=ground_field)),
        ("noon 7", Groebner.noonn(7, k=ground_field)),
        ("noon 8", Groebner.noonn(8, k=ground_field)),
        ("noon 9", Groebner.noonn(9, k=ground_field)),
        # ("noon 10", Groebner.noonn(10, k=ground_field)),
        # ("noon 11", Groebner.noonn(11, k=ground_field)),
        ("henrion 5", Groebner.henrion5(k=ground_field)),
        ("henrion 6", Groebner.henrion6(k=ground_field)),
        ("henrion 7", Groebner.henrion7(k=ground_field)),
        ("henrion 8", Groebner.henrion8(k=ground_field)),
        ("reimer 6", Groebner.reimern(6, k=ground_field)),
        ("reimer 7", Groebner.reimern(7, k=ground_field)),
        ("reimer 8", Groebner.reimern(8, k=ground_field))
        # ("reimer 9", Groebner.reimern(9, k=ground_field)),
        # ("chandra 11", Groebner.chandran(11, k=ground_field)),
        # ("chandra 12", Groebner.chandran(12, k=ground_field)),
        # ("chandra 13", Groebner.chandran(13, k=ground_field)),
        # ("chandra 14", Groebner.chandran(14, k=ground_field))
    ]

    (name="Integers modulo 2^30 + 3", field=ground_field, systems=systems)
end

function benchmark_set_2()
    ground_field = AbstractAlgebra.GF(1031)

    systems = [
        dummy_system("dummy", ground_field),
        ("cyclic 7", Groebner.cyclicn(7, k=ground_field)),
        ("cyclic 8", Groebner.cyclicn(8, k=ground_field)),
        ("cyclic 9", Groebner.cyclicn(9, k=ground_field)),
        # ("cyclic 10", Groebner.cyclicn(10, k=ground_field)),
        ("katsura 10", Groebner.katsuran(10, k=ground_field)),
        ("katsura 11", Groebner.katsuran(11, k=ground_field)),
        ("katsura 12", Groebner.katsuran(12, k=ground_field)),
        ("katsura 13", Groebner.katsuran(13, k=ground_field)),
        ("eco 11", Groebner.eco11(k=ground_field)),
        ("eco 12", Groebner.eco12(k=ground_field)),
        ("eco 13", Groebner.eco13(k=ground_field)),
        ("eco 14", Groebner.eco14(k=ground_field)),
        ("noon 7", Groebner.noonn(7, k=ground_field)),
        ("noon 8", Groebner.noonn(8, k=ground_field)),
        ("noon 9", Groebner.noonn(9, k=ground_field)),
        ("noon 10", Groebner.noonn(10, k=ground_field)),
        # ("noon 11", Groebner.noonn(11, k=ground_field)),
        ("henrion 5", Groebner.henrion5(k=ground_field)),
        ("henrion 6", Groebner.henrion6(k=ground_field)),
        ("henrion 7", Groebner.henrion7(k=ground_field)),
        ("henrion 8", Groebner.henrion8(k=ground_field)),
        ("reimer 6", Groebner.reimern(6, k=ground_field)),
        ("reimer 7", Groebner.reimern(7, k=ground_field)),
        ("reimer 8", Groebner.reimern(8, k=ground_field)),
        # ("reimer 9", Groebner.reimern(9, k=ground_field)),
        ("chandra 11", Groebner.chandran(11, k=ground_field)),
        ("chandra 12", Groebner.chandran(12, k=ground_field)),
        ("chandra 13", Groebner.chandran(13, k=ground_field)),
        ("chandra 14", Groebner.chandran(14, k=ground_field))
    ]

    (name="Integers modulo 1031", field=ground_field, systems=systems)
end

function benchmark_set_3()
    ground_field = AbstractAlgebra.QQ
    systems = [
        dummy_system("dummy", ground_field),
        ("cyclic 7", Groebner.cyclicn(7, k=ground_field)),
        ("cyclic 8", Groebner.cyclicn(8, k=ground_field)),
        # ("cyclic 9", Groebner.cyclicn(9, k=ground_field)),
        ("katsura 9", Groebner.katsuran(9, k=ground_field)),
        ("katsura 10", Groebner.katsuran(10, k=ground_field)),
        ("katsura 11", Groebner.katsuran(11, k=ground_field)),
        ("eco 10", Groebner.eco10(k=ground_field)),
        ("eco 11", Groebner.eco11(k=ground_field)),
        ("eco 12", Groebner.eco12(k=ground_field)),
        ("eco 13", Groebner.eco13(k=ground_field)),
        ("noon 7", Groebner.noonn(7, k=ground_field)),
        ("noon 8", Groebner.noonn(8, k=ground_field)),
        ("noon 9", Groebner.noonn(9, k=ground_field)),
        ("henrion 6", Groebner.henrion6(k=ground_field)),
        ("henrion 7", Groebner.henrion7(k=ground_field)),
        ("reimer 6", Groebner.reimern(6, k=ground_field)),
        ("reimer 7", Groebner.reimern(7, k=ground_field)),
        ("reimer 8", Groebner.reimern(8, k=ground_field)),
        ("chandra 9", Groebner.chandran(9, k=ground_field, internal_ordering=:degrevlex)),
        ("chandra 10", Groebner.chandran(10, k=ground_field, internal_ordering=:degrevlex)),
        ("chandra 11", Groebner.chandran(11, k=ground_field, internal_ordering=:degrevlex)),
        ("chandra 12", Groebner.chandran(12, k=ground_field, internal_ordering=:degrevlex)),
        ("chandra 13", Groebner.chandran(13, k=ground_field, internal_ordering=:degrevlex)),
        ("ipp", Groebner.ipp(k=ground_field, tol=0.0, internal_ordering=:degrevlex))
    ]

    (name="The rationals", field=ground_field, systems=systems)
end

function benchmark_set_4()
    ground_field = AbstractAlgebra.GF(2^30 + 3)
    systems = load_SIAN_all(ground=ground_field)

    (name="SIAN modulo 2^30 + 3", field=ground_field, systems=systems)
end

function benchmark_set_5()
    # TODO: this is not correct!!
    ground_field = AbstractAlgebra.GF(2^30 + 3)
    names = [
        "mq_n10_m20_p2_s0",
        "mq_n10_m20_p2_s1",
        "mq_n10_m7_p2_s0",
        "mq_n15_m10_p2_s0",
        "mq_n15_m30_p2_s0",
        "mq_n24_m16_p31_s0",
        "mq_n34_m68_p31_s0"
    ]
    systems = [(name, load_MQ_problem(name)) for name in names]

    (name="MQ", field=ground_field, systems=systems)
end

function benchmark_set_6()
    ground_field = AbstractAlgebra.QQ
    names = ["SIWR", "SEAIJRC"]
    systems = [(name, load_SI_problem(name)) for name in names]

    (name="SI", field=ground_field, systems=systems)
end

function benchmark_set_7()
    ground_field = AbstractAlgebra.QQ
    systems = [
        ("chandra 2", Groebner.chandran(2, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 3", Groebner.chandran(3, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 4", Groebner.chandran(4, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 5", Groebner.chandran(5, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 6", Groebner.chandran(6, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 7", Groebner.chandran(7, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 8", Groebner.chandran(8, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 9", Groebner.chandran(9, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 10", Groebner.chandran(10, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 11", Groebner.chandran(11, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 12", Groebner.chandran(12, internal_ordering=:degrevlex, k=ground_field)),
        ("boon", Groebner.boon(internal_ordering=:degrevlex, k=ground_field)),
        ("rps10", Groebner.rps10(internal_ordering=:degrevlex, k=ground_field)),
        ("ipp", Groebner.ipp(internal_ordering=:degrevlex, k=ground_field))
    ]

    (name="HC", field=ground_field, systems=systems)
end

function benchmark_set_8()
    ground_field = AbstractAlgebra.GF(2^30 + 3)
    systems = [
        ("chandra 5", Groebner.chandran(5, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 6", Groebner.chandran(6, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 7", Groebner.chandran(7, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 8", Groebner.chandran(8, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 9", Groebner.chandran(9, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 10", Groebner.chandran(10, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 11", Groebner.chandran(11, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 12", Groebner.chandran(12, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 13", Groebner.chandran(13, internal_ordering=:degrevlex, k=ground_field)),
        ("chandra 14", Groebner.chandran(14, internal_ordering=:degrevlex, k=ground_field))
    ]

    (name="HC modulo 2^30 + 3", field=ground_field, systems=systems)
end

function benchmark_set_9()
    ground_field = AbstractAlgebra.QQ
    systems =
        [("BIOMD0000000103", BIOMD0000000103()), ("BIOMD0000000123", BIOMD0000000123())]
    (name="BIOMD, the rationals", field=ground_field, systems=systems)
end

function benchmark_set_10()
    ground_field = AbstractAlgebra.QQ
    systems = load_SIAN_all(ground=ground_field)

    (name="SIAN, the rationals", field=ground_field, systems=systems)
end

function benchmark_set_11()
    ground_field = AbstractAlgebra.QQ
    systems = [
    # ("hexapod", Groebner.hexapod(k=ground_field)), 
        ("ipp", Groebner.ipp(k=ground_field))
    # ("alea6", Groebner.alea6(k=ground_field))
    ]

    (name="ipp, the rationals", field=ground_field, systems=systems)
end

function benchmark_set_12()
    ground_field = AbstractAlgebra.GF(2^30 + 3)
    systems = [
        dummy_system("dummy", ground_field),
        # ("cyclic 9", Groebner.cyclicn(9, k=ground_field)),
        # ("katsura 12", Groebner.katsuran(12, k=ground_field)),
        # ("eco 14", Groebner.eco14(k=ground_field)),
        # ("noon 10", Groebner.noonn(10, k=ground_field)),
        ("reimer 9", Groebner.reimern(9, k=ground_field)),
        ("eco 15", Groebner.eco15(k=ground_field))
    ]

    (name="2^30+3, larger", field=ground_field, systems=systems)
end

function benchmark_set_13()
    ground_field = GF(2^30 + 3)
    systems = [
        dummy_system("dummy", ground_field),
        ("ipp", Groebner.hexapod(internal_ordering=:degrevlex, k=ground_field)),
        ("ipp", Groebner.ipp(internal_ordering=:degrevlex, k=ground_field)),
        # ("SIWR", load_SI_problem("SIWR", k=ground_field)),
        # ("SEAIJRC", load_SI_problem("SEAIJRC", k=ground_field)),
        ("BIOMD0000000103", BIOMD0000000103(k=ground_field)),
        ("BIOMD0000000123", BIOMD0000000123(k=ground_field)),
        ("Cholera", Cholera(k=ground_field)),
        ("HIV2", HIV2(k=ground_field)),
        ("NFkB", NFkB(k=ground_field)),
        ("Pharm_with_weights", Pharm_with_weights(k=ground_field)),
        ("SEIRP", SEIRP(k=ground_field))
    ]

    (name="Learn-apply, other", field=ground_field, systems=systems)
end

function benchmark_set_14()
    ground_field = GF(2^30 + 3)
    systems = [
        dummy_system("dummy", ground_field),
        ("yang1", Groebner.yang1(internal_ordering=:degrevlex, k=ground_field)),
        ("bayes148", Groebner.bayes148(internal_ordering=:degrevlex, k=ground_field)),
        ("gametwo2", Groebner.gametwo2(internal_ordering=:degrevlex, k=ground_field)),
        ("jason210", Groebner.jason210(internal_ordering=:degrevlex, k=ground_field)),
        ("alea6", Groebner.alea6(internal_ordering=:degrevlex, k=ground_field)),
        ("mayr42", Groebner.mayr42(internal_ordering=:degrevlex, k=ground_field))
    ]

    (name="Other, modulo 2^30 + 3", field=ground_field, systems=systems)
end

function benchmark_set_15()
    ground_field = GF(2^30 + 3)
    systems = [
        dummy_system("dummy", ground_field),
        ("Cholera", Cholera(k=ground_field)),
        ("HIV2", HIV2(k=ground_field)),
        ("NFkB (w.)", NFkB_with_weights(k=ground_field)),
        # ("Pharm (w.)", Pharm_with_weights(k=ground_field)),
        ("Goodwin (w.)", Goodwin_with_weights(k=ground_field))
        # ("SEIRP", SEIRP(k=ground_field)),
    ]

    (name="SIAN, 2^30+3", field=ground_field, systems=systems)
end

function benchmark_set_16()
    ground_field = AbstractAlgebra.GF(2^30 + 3)
    systems = [
        dummy_system("dummy", ground_field),
        # ("cyclic 9", Groebner.cyclicn(9, k=ground_field)),
        # ("katsura 12", Groebner.katsuran(12, k=ground_field)),
        # ("eco 14", Groebner.eco14(k=ground_field)),
        ("noon 10", Groebner.noonn(10, k=ground_field))
        # ("reimer 9", Groebner.reimern(9, k=ground_field)),
        # ("eco 15", Groebner.eco15(k=ground_field)),
    ]

    (name="2^30+3, larger x2", field=ground_field, systems=systems)
end

function benchmark_set_17()
    ground_field = QQ
    systems = [
        dummy_system("dummy", ground_field),
        ("Cholera", Cholera(k=ground_field)),
        ("HIV2", HIV2(k=ground_field)),
        ("NFkB (w.)", NFkB_with_weights(k=ground_field)),
        ("Pharm (w.)", Pharm_with_weights(k=ground_field)),
        ("Goodwin (w.)", Goodwin_with_weights(k=ground_field)),
        ("SEIRP", SEIRP(k=ground_field))
    ]

    (name="SIAN, QQ", field=ground_field, systems=systems)
end

function benchmark_set_18()
    ground_field = AbstractAlgebra.GF(2^30 + 3)

    systems = [
        dummy_system("dummy", ground_field),
        ("cyclic 7", Groebner.cyclicn(7, k=ground_field)),
        ("cyclic 8", Groebner.cyclicn(8, k=ground_field)),
        ("cyclic 9", Groebner.cyclicn(9, k=ground_field)),
        # ("cyclic 10", Groebner.cyclicn(10, k=ground_field)),
        ("katsura 10", Groebner.katsuran(10, k=ground_field)),
        ("katsura 11", Groebner.katsuran(11, k=ground_field)),
        ("katsura 12", Groebner.katsuran(12, k=ground_field)),
        ("katsura 13", Groebner.katsuran(13, k=ground_field)),
        # ("katsura 14", Groebner.katsuran(14, k=ground_field)),
        ("eco 11", Groebner.eco11(k=ground_field)),
        ("eco 12", Groebner.eco12(k=ground_field)),
        ("eco 13", Groebner.eco13(k=ground_field)),
        ("eco 14", Groebner.eco14(k=ground_field)),
        ("noon 7", Groebner.noonn(7, k=ground_field)),
        ("noon 8", Groebner.noonn(8, k=ground_field)),
        ("noon 9", Groebner.noonn(9, k=ground_field)),
        # ("noon 10", Groebner.noonn(10, k=ground_field)),
        # ("noon 11", Groebner.noonn(11, k=ground_field)),
        ("henrion 5", Groebner.henrion5(k=ground_field)),
        ("henrion 6", Groebner.henrion6(k=ground_field)),
        ("henrion 7", Groebner.henrion7(k=ground_field)),
        ("henrion 8", Groebner.henrion8(k=ground_field)),
        ("reimer 6", Groebner.reimern(6, k=ground_field)),
        ("reimer 7", Groebner.reimern(7, k=ground_field)),
        ("reimer 8", Groebner.reimern(8, k=ground_field))
        # ("reimer 9", Groebner.reimern(9, k=ground_field)),
        # ("chandra 11", Groebner.chandran(11, k=ground_field)),
        # ("chandra 12", Groebner.chandran(12, k=ground_field)),
        # ("chandra 13", Groebner.chandran(13, k=ground_field)),
        # ("chandra 14", Groebner.chandran(14, k=ground_field))
    ]

    (name="Integers modulo 2^30 + 3", field=ground_field, systems=systems)
end
