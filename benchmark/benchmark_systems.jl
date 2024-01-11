import AbstractAlgebra, Groebner

function get_benchmark_suite(id)
    if id == 0
        benchmark_set_0()
    elseif id == 1
        benchmark_set_1()
    elseif id == 2
        benchmark_set_2()
    elseif id == 3
        benchmark_set_3()
    elseif id == 4
        benchmark_set_4()
    end
end

function dummy_system(name, ground_field)
    ring, (x, y) = AbstractAlgebra.polynomial_ring(ground_field, ["x", "y"])
    (name, [x^2 + y^2 - 1, x + y])
end

function benchmark_set_0()
    (
        name="dummy benchmark set",
        field=AbstractAlgebra.QQ,
        systems=[
            dummy_system("dummy 1", AbstractAlgebra.QQ),
            dummy_system("dummy 2", AbstractAlgebra.QQ)
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
        ("cyclic 10", Groebner.cyclicn(10, k=ground_field)),
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
        ("henrion 5", Groebner.henrion5(k=ground_field)),
        ("henrion 6", Groebner.henrion6(k=ground_field)),
        ("henrion 7", Groebner.henrion7(k=ground_field)),
        ("henrion 8", Groebner.henrion8(k=ground_field)),
        ("reimer 6", Groebner.reimern(6, k=ground_field)),
        ("reimer 7", Groebner.reimern(7, k=ground_field)),
        ("reimer 8", Groebner.reimern(8, k=ground_field)),
        ("reimer 9", Groebner.reimern(9, k=ground_field))
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
        ("cyclic 10", Groebner.cyclicn(10, k=ground_field)),
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
        ("henrion 5", Groebner.henrion5(k=ground_field)),
        ("henrion 6", Groebner.henrion6(k=ground_field)),
        ("henrion 7", Groebner.henrion7(k=ground_field)),
        ("reimer 6", Groebner.reimern(6, k=ground_field)),
        ("reimer 7", Groebner.reimern(7, k=ground_field)),
        ("reimer 8", Groebner.reimern(8, k=ground_field)),
        ("reimer 9", Groebner.reimern(9, k=ground_field))
    ]

    (name="Integers modulo 1031", field=ground_field, systems=systems)
end

function benchmark_set_3()
    ground_field = AbstractAlgebra.QQ
    systems = [
        dummy_system("dummy", ground_field),
        ("cyclic 7", Groebner.cyclicn(7, k=ground_field)),
        ("cyclic 8", Groebner.cyclicn(8, k=ground_field)),
        ("katsura 9", Groebner.katsuran(9, k=ground_field)),
        ("katsura 10", Groebner.katsuran(10, k=ground_field)),
        ("katsura 11", Groebner.katsuran(11, k=ground_field)),
        ("eco 10", Groebner.eco10(k=ground_field)),
        ("eco 11", Groebner.eco11(k=ground_field)),
        ("eco 12", Groebner.eco12(k=ground_field)),
        ("noon 8", Groebner.noonn(8, k=ground_field)),
        ("noon 9", Groebner.noonn(9, k=ground_field)),
        ("henrion 6", Groebner.henrion6(k=ground_field)),
        ("henrion 7", Groebner.henrion7(k=ground_field)),
        ("reimer 6", Groebner.reimern(6, k=ground_field)),
        ("reimer 7", Groebner.reimern(7, k=ground_field)),
        ("reimer 8", Groebner.reimern(8, k=ground_field))
    ]

    (name="The rationals", field=ground_field, systems=systems)
end

function benchmark_set_4()
    ground_field = AbstractAlgebra.GF(2^30 + 3)
    systems = [
        dummy_system("dummy", ground_field),
        ("NF-kB", NF_kB(k=ground_field)),
        ("NF-kB, weights", NF_kB_with_weights(k=ground_field)),
        ("Chol", chol(k=ground_field)),
        ("Chol, 1 out", chol_1_out(k=ground_field)),
        ("Chol, weights", chol_with_weights(k=ground_field)),
        ("Pharm", Pharm(k=ground_field)),
        ("Pharm_with_weights", Pharm(k=ground_field))
    ]

    (name="SIAN", field=ground_field, systems=systems)
end
