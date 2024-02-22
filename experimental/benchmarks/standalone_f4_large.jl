
# using Groebner

import AbstractAlgebra
using BenchmarkTools
using Logging
# global_logger(ConsoleLogger(stderr, Logging.Error))

Groebner.logging_enabled() = false
Groebner.invariants_enabled() = false

function benchmark_system_my(system)
    # system = Groebner.change_ordering(system, :degrevlex)

    Groebner.groebner(system)

    t = Inf
    for i in 1:5
        timing = @timed Groebner.groebner(system)
        t = min(t, timing.time)
    end
    t
end

function run_f4_ff_degrevlex_benchmarks(ground)
    systems = [
        ("cyclic 7", Groebner.cyclicn(7, k=ground, internal_ordering=:degrevlex)),
        ("cyclic 8", Groebner.cyclicn(8, k=ground, internal_ordering=:degrevlex)),
        # ("cyclic 9", Groebner.cyclicn(9, k=ground, internal_ordering=:degrevlex)),
        ("katsura 9", Groebner.katsuran(9, k=ground, internal_ordering=:degrevlex)),
        ("katsura 10", Groebner.katsuran(10, k=ground, internal_ordering=:degrevlex)),
        ("katsura 11", Groebner.katsuran(11, k=ground, internal_ordering=:degrevlex)),
        # ("katsura 12", Groebner.katsuran(12, k=ground, internal_ordering=:degrevlex)),
        ("noon 7", Groebner.noonn(7, k=ground, internal_ordering=:degrevlex)),
        ("noon 8", Groebner.noonn(8, k=ground, internal_ordering=:degrevlex)),
        ("noon 9", Groebner.noonn(9, k=ground, internal_ordering=:degrevlex)),
        ("eco 11", Groebner.eco11(k=ground, internal_ordering=:degrevlex)),
        ("eco 12", Groebner.eco12(k=ground, internal_ordering=:degrevlex)),
        ("eco 13", Groebner.eco13(k=ground, internal_ordering=:degrevlex))
    ]

    for (name, system) in systems
        print("$name\t--\t")
        t = benchmark_system_my(system)
        println("$t s")
    end
end

ground = AbstractAlgebra.GF(2^30 + 3)
run_f4_ff_degrevlex_benchmarks(ground)

#=
inline      : -,
assume Zp   : no
assume empty: no 
cyclic 7        --      0.0965658 s
cyclic 8        --      1.4198005 s
katsura 9       --      0.1470269 s
katsura 10      --      0.8665618 s
katsura 11      --      5.2572539 s
noon 7  --      0.1516447 s
noon 8  --      1.2835755 s
noon 9  --      12.139458 s
eco 11  --      0.3524228 s
eco 12  --      2.1887455 s
eco 13  --      9.0118321 s

inline      : -,
assume Zp   : yes
assume empty: no 
cyclic 7        --      0.0894644 s
cyclic 8        --      1.2891848 s
katsura 9       --      0.1435553 s
katsura 10      --      0.8241746 s
katsura 11      --      4.9737961 s
noon 7  --      0.1593515 s
noon 8  --      1.3084553 s
noon 9  --      12.112527 s
eco 11  --      0.3444795 s
eco 12  --      2.1214438 s
eco 13  --      8.6668539 s

inline      : -,
assume Zp   : yes
assume empty: yes
cyclic 7        --      0.0875634 s
cyclic 8        --      1.29229 s
katsura 9       --      0.1403045 s
katsura 10      --      0.8176833 s
katsura 11      --      4.9808602 s
noon 7  --      0.152636 s
noon 8  --      1.301324 s
noon 9  --      11.99392 s
eco 11  --      0.3462996 s
eco 12  --      2.091162 s
eco 13  --      8.5554278 s

inline      : yes,
assume Zp   : yes
assume empty: yes
cyclic 7        --      0.0892941 s
cyclic 8        --      1.3551499 s
katsura 9       --      0.1469254 s
katsura 10      --      0.8492118 s
katsura 11      --      5.1888122 s
noon 7  --      0.1558268 s
noon 8  --      1.3372536 s
noon 9  --      12.066885 s
eco 11  --      0.345839 s
eco 12  --      2.1880031 s
eco 13  --      8.9609859 s
=#
