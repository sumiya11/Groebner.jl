
if !isdefined(Main, :Groebner)
    import Groebner
end

import AbstractAlgebra
using BenchmarkTools
using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

function my_system_coeffs_1(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    gb = Groebner.groebner(system, reduced=false)

    aver = [(0, 0) for i in 1:length(gb)]
    for (i, poly) in enumerate(gb)
        num = sum([Base.GMP.MPZ.sizeinbase(numerator(c), 2) for c in AbstractAlgebra.coefficients(poly)])
        den = sum([Base.GMP.MPZ.sizeinbase(denominator(c), 2) for c in AbstractAlgebra.coefficients(poly)])

        aver[i] = (num, den)
    end

    Groebner.printall()

    # println(aver)

    println("Numerators: ", sum(x[1] for x in aver))
    println("Denominators: ", sum(x[2] for x in aver))

end

function my_system_coeffs_2(system)
    system = Groebner.change_ordering(system, :degrevlex)
    Groebner.groebner([system[1]])

    # gb = Groebner.groebner(system)
    # println("length = $(length(gb))")
    # println("degree = $(maximum(AbstractAlgebra.total_degree, gb))")

    gb = Groebner.groebner(system, reduced=false)

    aver = [(0, 0) for i in 1:length(gb)]
    for (i, poly) in enumerate(gb)
        @warn "new poly" i
        println([(Base.GMP.MPZ.sizeinbase(numerator(c), 2), Base.GMP.MPZ.sizeinbase(denominator(c), 2)) for c in AbstractAlgebra.coefficients(poly)])
    end

    Groebner.printall()

    # println(aver)

    println("Numerators: ", sum(x[1] for x in aver))
    println("Denominators: ", sum(x[2] for x in aver))

end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        #("cyclic 7", Groebner.cyclicn(7, ground=ground)),
        #("noon 8"    ,Groebner.noonn(8, ground=ground)),
        # ("eco 11"    ,Groebner.eco11(ground=ground)),
        ("katsura 8"    ,Groebner.katsuran(8, ground=ground)),
        #("root 11"    ,Groebner.rootn(11, ground=ground)),
        ]

    println()
    for (name, system) in systems
        println("-"^20)
        println("$name")
        # benchmark_system_my(system)
        my_system_coeffs_2(system)
        println("-"^20)
    end
end

ground = AbstractAlgebra.QQ
run_f4_ff_degrevlex_benchmarks(ground)

#=
cyclic 11
  50.711 ms (191971 allocations: 23.54 MiB)
noon 7
  910.447 ms (2087052 allocations: 298.79 MiB)
eco 10
  1.484 s (1994909 allocations: 417.28 MiB)
katsura 8
  4.509 s (2982135 allocations: 499.90 MiB)
=#
