# This script is not to be run directly. See runtests.jl.

stopwatch = time_ns()

# TTFX
t1 = @timed using Groebner
using AbstractAlgebra

R, (x, y) = polynomial_ring(QQ, ["x", "y"], internal_ordering=:degrevlex)
sys = [x^2 + y, x * y^2]
t2 = @timed groebner(sys)

R, (x, y) = polynomial_ring(GF(2^30 + 3), ["x", "y"])
sys = [x^2 + y, x * y^2]
t3 = @timed groebner(sys)

suite = []

push!(suite, (problem_name="using Groebner", type=:time, result=[t1.time]))
push!(suite, (problem_name="groebner,ttfx,qq drl", type=:time, result=[t2.time]))
push!(suite, (problem_name="groebner,ttfx,zp lex", type=:time, result=[t3.time]))

# Allocations
k = Groebner.Examples.katsuran(7, internal_ordering=:degrevlex)
groebner(k)
a1 = @allocated groebner(k)

k = Groebner.Examples.katsuran(7, k=GF(2^30 + 3), internal_ordering=:degrevlex)
groebner(k)
a2 = @allocated groebner(k)

push!(suite, (problem_name="AA,bytes,QQ,katsura-7", type=:allocs, result=[a1]))
push!(suite, (problem_name="AA,bytes,2^30+3,katsura-7", type=:allocs, result=[a2]))

k = Groebner.Examples.katsuran(7, internal_ordering=:degrevlex)
a1 = @allocations groebner(k)

k = Groebner.Examples.katsuran(7, k=GF(2^30 + 3), internal_ordering=:degrevlex)
a2 = @allocations groebner(k)

push!(suite, (problem_name="AA,allocs,QQ,katsura-7", type=:allocs, result=[a1]))
push!(suite, (problem_name="AA,allocs,2^30+3,katsura-7", type=:allocs, result=[a2]))

k = Groebner.Examples.cyclicn(8, k=GF(2^30 + 3), internal_ordering=:degrevlex)
groebner(k)
a3 = @allocations groebner(k)
push!(suite, (problem_name="AA,allocs,2^30+3,cyclic-8", type=:allocs, result=[a3]))

k = Groebner.Examples.katsuran(11, k=GF(2^30 + 3), internal_ordering=:degrevlex)
groebner(k)
a4 = @allocations groebner(k)
push!(suite, (problem_name="AA,allocs,2^30+3,katsura-11", type=:allocs, result=[a4]))

# Runtime
import Primes
import Nemo

function nemo_make_prime_finite_field(p)
    if p < typemax(UInt)
        Nemo.fpField(convert(UInt, p), false)
    else
        Nemo.FpField(Nemo.ZZRingElem(p), false)
    end
end

function compute_gb(system, trials=7; kws...)
    times = []
    for _ in 1:trials
        GC.gc()
        time = @elapsed groebner(system; kws...)
        push!(times, time)
    end
    times
end

# Compute Groebner bases over integers modulo a large prime
push!(
    suite,
    (
        problem_name="groebner,AA,2^31-1,katsura 9",
        type=:time,
        result=compute_gb(
            Groebner.Examples.katsuran(9, internal_ordering=:degrevlex, k=GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,2^31-1,katsura 10",
        type=:time,
        result=compute_gb(
            Groebner.Examples.katsuran(10, internal_ordering=:degrevlex, k=GF(2^31 - 1)),
            5
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,2^31-1,katsura 11",
        type=:time,
        result=compute_gb(
            Groebner.Examples.katsuran(11, internal_ordering=:degrevlex, k=GF(2^31 - 1)),
            3
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,117,yang1",
        type=:time,
        result=compute_gb(
            Groebner.Examples.yang1(internal_ordering=:degrevlex, k=GF(117)),
            2
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,2^27+29,katsura 10",
        type=:time,
        result=compute_gb(
            Groebner.Examples.katsuran(10, internal_ordering=:degrevlex, k=GF(2^27 + 29)),
            5
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,2^27+29,cyclic 8",
        type=:time,
        result=compute_gb(
            Groebner.Examples.cyclicn(8, internal_ordering=:degrevlex, k=GF(2^27 + 29)),
            5
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,2^31-1,cyclic 8",
        type=:time,
        result=compute_gb(
            Groebner.Examples.cyclicn(8, internal_ordering=:degrevlex, k=GF(2^31 - 1)),
            5
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner,Nemo,2^31-1,cyclic 8",
        type=:time,
        result=compute_gb(
            Groebner.Examples.cyclicn(
                8,
                internal_ordering=:degrevlex,
                k=nemo_make_prime_finite_field(2^31 - 1)
            ),
            5
        )
    )
)

function n_variable_set(n; internal_ordering=:degrevlex, k=GF(2^31 - 1))
    R, x = polynomial_ring(k, ["x$i" for i in 1:n], internal_ordering=internal_ordering)
    f = [sum(prod(x[i:(n - kk)], init=1) for i in 1:(kk + 1)) for kk in 0:(n - 1)]
    f
end

push!(
    suite,
    (
        problem_name="groebner,AA,100 vars drl",
        type=:time,
        result=compute_gb(
            n_variable_set(100, internal_ordering=:degrevlex, k=GF(2^31 - 1))
        )
    )
)

push!(
    suite,
    (
        problem_name="groebner,threaded,AA,2^31-1,cyclic 8",
        type=:time,
        result=compute_gb(
            Groebner.Examples.cyclicn(8, internal_ordering=:degrevlex, k=GF(2^31 - 1)),
            5,
            threaded=:yes
        )
    )
)

function learn_and_apply(system)
    times = []
    trials = 7
    graph, gb = groebner_learn(system)
    for _ in 1:trials
        GC.gc()
        time2 = @elapsed groebner_apply!(graph, system)
        push!(times, time2)
    end
    times
end

# Learn and apply
push!(
    suite,
    (
        problem_name="apply,AA,2^31-1,cyclic 7",
        type=:time,
        result=learn_and_apply(
            Groebner.Examples.cyclicn(7, internal_ordering=:degrevlex, k=GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="apply,Nemo,2^31-1,cyclic 7",
        type=:time,
        result=learn_and_apply(
            Groebner.Examples.cyclicn(
                7,
                internal_ordering=:degrevlex,
                k=nemo_make_prime_finite_field(2^31 - 1)
            )
        )
    )
)
push!(
    suite,
    (
        problem_name="apply,AA,2^31-1,katsura 10",
        type=:time,
        result=learn_and_apply(
            Groebner.Examples.katsuran(10, internal_ordering=:degrevlex, k=GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="apply,AA,2^27+29,katsura 10",
        type=:time,
        result=learn_and_apply(
            Groebner.Examples.katsuran(10, internal_ordering=:degrevlex, k=GF(2^27 + 29))
        )
    )
)
push!(
    suite,
    (
        problem_name="apply,Nemo,2^31-1,katsura 10",
        type=:time,
        result=learn_and_apply(
            Groebner.Examples.katsuran(
                10,
                internal_ordering=:degrevlex,
                k=nemo_make_prime_finite_field(2^31 - 1)
            )
        )
    )
)

# Compute Groebner bases over the rationals
push!(
    suite,
    (
        problem_name="groebner,AA,QQ,katsura 8",
        type=:time,
        result=compute_gb(
            Groebner.Examples.katsuran(8, internal_ordering=:degrevlex, k=QQ)
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner,Nemo,QQ,katsura 8",
        type=:time,
        result=compute_gb(
            Groebner.Examples.katsuran(8, internal_ordering=:degrevlex, k=Nemo.QQ)
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,QQ,eco 10",
        type=:time,
        result=compute_gb(Groebner.Examples.eco10(internal_ordering=:degrevlex, k=QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,QQ,cyclic 7",
        type=:time,
        result=compute_gb(Groebner.Examples.cyclicn(7, internal_ordering=:degrevlex, k=QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,QQ,noon 8",
        type=:time,
        result=compute_gb(
            Groebner.Examples.noonn(8, internal_ordering=:degrevlex, k=QQ),
            3
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,QQ,hexapod",
        type=:time,
        result=compute_gb(Groebner.Examples.hexapod(internal_ordering=:degrevlex, k=QQ), 3)
    )
)

function multimodular_gb_problem(nbits; np=AbstractAlgebra)
    R, (x1, x2, x3, x4) =
        polynomial_ring(np.QQ, ["x1", "x2", "x3", "x4"], internal_ordering=:degrevlex)
    nbits_per_prime = 31
    nprimes = max(div(nbits, nbits_per_prime), 1)
    N = prod(map(BigInt, Primes.nextprimes(2^31 - 100, nprimes)))
    system = [
        x1 + x2 + x3 + x4,
        x1 * x2 + x1 * x3 + x1 * x4 + x2 * x3 + x2 * x4 + x3 * x4,
        x1 * x2 * x3 + x1 * x2 * x4 + x1 * x3 * x4 + x2 * x3 * x4,
        x1 * x2 * x3 * x4 + N
    ]
    system
end

push!(
    suite,
    (
        problem_name="groebner,AA,QQ,1000-bit output",
        type=:time,
        result=compute_gb(multimodular_gb_problem(1000, np=AbstractAlgebra))
    )
)
push!(
    suite,
    (
        problem_name="groebner,Nemo,QQ,1000-bit output",
        type=:time,
        result=compute_gb(multimodular_gb_problem(1000, np=Nemo))
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,QQ,10_000-bit output",
        type=:time,
        result=compute_gb(multimodular_gb_problem(10000, np=AbstractAlgebra))
    )
)
push!(
    suite,
    (
        problem_name="groebner,Nemo,QQ,10_000-bit output",
        type=:time,
        result=compute_gb(multimodular_gb_problem(10000, np=Nemo))
    )
)
push!(
    suite,
    (
        problem_name="groebner,AA,QQ,100_000-bit output",
        type=:time,
        result=compute_gb(multimodular_gb_problem(100000, np=AbstractAlgebra), 3)
    )
)
push!(
    suite,
    (
        problem_name="groebner,Nemo,QQ,100_000-bit output",
        type=:time,
        result=compute_gb(multimodular_gb_problem(100000, np=Nemo), 3)
    )
)

function compute_normalforms(system)
    R = parent(system[1])
    gb = Groebner.groebner(system)
    times = []
    trials = 7
    for _ in 1:trials
        GC.gc()
        time = @elapsed begin
            n1 = normalform(gb, system)
            n2 = normalform(gb, gb)
        end
        push!(times, time)
    end
    times
end

# Normal forms
push!(
    suite,
    (
        problem_name="normalform,AA,103,cyclic 8",
        type=:time,
        result=compute_normalforms(
            Groebner.Examples.cyclicn(8, internal_ordering=:degrevlex, k=GF(103))
        )
    )
)
push!(
    suite,
    (
        problem_name="normalform,Nemo,103,cyclic 8",
        type=:time,
        result=compute_normalforms(
            Groebner.Examples.cyclicn(
                8,
                internal_ordering=:degrevlex,
                k=nemo_make_prime_finite_field(103)
            )
        )
    )
)
push!(
    suite,
    (
        problem_name="normalform,AA,QQ,katsura 9",
        type=:time,
        result=compute_normalforms(
            Groebner.Examples.katsuran(9, internal_ordering=:degrevlex, k=QQ)
        )
    )
)

stopwatch = time_ns() - stopwatch
push!(suite, (problem_name="total", type=:time, result=[stopwatch / 1e9]))
