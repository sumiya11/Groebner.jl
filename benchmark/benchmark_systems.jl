import AbstractAlgebra, Groebner

function get_benchmark(code)
    if code == 0
        benchmark_set_0()
    elseif code == 1
        benchmark_set_1()
    elseif code == 2
        benchmark_set_2()
    elseif code == 3
        benchmark_set_3()
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
        ("cyclic 7", Groebner.cyclicn(7, ground=ground_field)),
        ("cyclic 8", Groebner.cyclicn(8, ground=ground_field)),
        ("cyclic 9", Groebner.cyclicn(9, ground=ground_field)),
        ("cyclic 10", Groebner.cyclicn(10, ground=ground_field)),
        ("katsura 10", Groebner.katsuran(10, ground=ground_field)),
        ("katsura 11", Groebner.katsuran(11, ground=ground_field)),
        ("katsura 12", Groebner.katsuran(12, ground=ground_field)),
        ("katsura 13", Groebner.katsuran(13, ground=ground_field)),
        ("eco 11", Groebner.eco11(ground=ground_field)),
        ("eco 12", Groebner.eco12(ground=ground_field)),
        ("eco 13", Groebner.eco13(ground=ground_field)),
        ("eco 14", Groebner.eco14(ground=ground_field)),
        ("noon 7", Groebner.noonn(7, ground=ground_field)),
        ("noon 8", Groebner.noonn(8, ground=ground_field)),
        ("noon 9", Groebner.noonn(9, ground=ground_field)),
        ("noon 10", Groebner.noonn(10, ground=ground_field)),
        ("henrion 5", Groebner.henrion5(ground=ground_field)),
        ("henrion 6", Groebner.henrion6(ground=ground_field)),
        ("henrion 7", Groebner.henrion7(ground=ground_field)),
        ("reimer 6", Groebner.reimern(6, ground=ground_field)),
        ("reimer 7", Groebner.reimern(7, ground=ground_field)),
        ("reimer 8", Groebner.reimern(8, ground=ground_field)),
        ("reimer 9", Groebner.reimern(9, ground=ground_field))
    ]

    (name="Integers modulo 2^30 + 3", field=ground_field, systems=systems)
end

function benchmark_set_2()
    ground_field = AbstractAlgebra.GF(1031)

    systems = [
        dummy_system("dummy", ground_field),
        ("cyclic 7", Groebner.cyclicn(7, ground=ground_field)),
        ("cyclic 8", Groebner.cyclicn(8, ground=ground_field)),
        ("cyclic 9", Groebner.cyclicn(9, ground=ground_field)),
        ("cyclic 10", Groebner.cyclicn(10, ground=ground_field)),
        ("katsura 10", Groebner.katsuran(10, ground=ground_field)),
        ("katsura 11", Groebner.katsuran(11, ground=ground_field)),
        ("katsura 12", Groebner.katsuran(12, ground=ground_field)),
        ("katsura 13", Groebner.katsuran(13, ground=ground_field)),
        ("eco 11", Groebner.eco11(ground=ground_field)),
        ("eco 12", Groebner.eco12(ground=ground_field)),
        ("eco 13", Groebner.eco13(ground=ground_field)),
        ("eco 14", Groebner.eco14(ground=ground_field)),
        ("noon 7", Groebner.noonn(7, ground=ground_field)),
        ("noon 8", Groebner.noonn(8, ground=ground_field)),
        ("noon 9", Groebner.noonn(9, ground=ground_field)),
        ("noon 10", Groebner.noonn(10, ground=ground_field)),
        ("henrion 5", Groebner.henrion5(ground=ground_field)),
        ("henrion 6", Groebner.henrion6(ground=ground_field)),
        ("henrion 7", Groebner.henrion7(ground=ground_field)),
        ("reimer 6", Groebner.reimern(6, ground=ground_field)),
        ("reimer 7", Groebner.reimern(7, ground=ground_field)),
        ("reimer 8", Groebner.reimern(8, ground=ground_field)),
        ("reimer 9", Groebner.reimern(9, ground=ground_field))
    ]

    (name="Integers modulo 1031", field=ground_field, systems=systems)
end

function benchmark_set_3()
    ground_field = AbstractAlgebra.QQ
    systems = [
        dummy_system("dummy", ground_field),
        ("cyclic 7", Groebner.cyclicn(7, ground=ground_field)),
        ("cyclic 8", Groebner.cyclicn(8, ground=ground_field)),
        ("katsura 9", Groebner.katsuran(9, ground=ground_field)),
        ("katsura 10", Groebner.katsuran(10, ground=ground_field)),
        ("katsura 11", Groebner.katsuran(11, ground=ground_field)),
        ("eco 10", Groebner.eco10(ground=ground_field)),
        ("eco 11", Groebner.eco11(ground=ground_field)),
        ("eco 12", Groebner.eco12(ground=ground_field)),
        ("noon 8", Groebner.noonn(8, ground=ground_field)),
        ("noon 9", Groebner.noonn(9, ground=ground_field)),
        ("henrion 6", Groebner.henrion6(ground=ground_field)),
        ("henrion 7", Groebner.henrion7(ground=ground_field)),
        ("reimer 6", Groebner.reimern(6, ground=ground_field)),
        ("reimer 7", Groebner.reimern(7, ground=ground_field)),
        ("reimer 8", Groebner.reimern(8, ground=ground_field))
    ]

    (name="The rationals", field=ground_field, systems=systems)
end
