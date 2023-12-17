# This script is not intended to be run directly.
# Instead, see runtests.jl.

using Groebner
using AbstractAlgebra
import Primes
import Nemo

suite = []

function compute_gb(system, trials=7)
    times = []
    for _ in 1:trials
        GC.gc()
        time = @elapsed groebner(system)
        push!(times, time)
    end
    minimum(times)
end

# Compute Groebner bases over integers modulo a large prime
problem = (
    problem_name="groebner, AA, GF(2^31-1), katsura 5",
    result=compute_gb(Groebner.katsuran(5, ordering=:degrevlex, ground=GF(2^31 - 1)))
)
push!(suite, problem)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), katsura 5",
        result=compute_gb(Groebner.katsuran(5, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), katsura 6",
        result=compute_gb(Groebner.katsuran(6, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), katsura 8",
        result=compute_gb(Groebner.katsuran(8, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), katsura 10",
        result=compute_gb(
            Groebner.katsuran(10, ordering=:degrevlex, ground=GF(2^31 - 1)),
            5
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^27+29), katsura 10",
        result=compute_gb(
            Groebner.katsuran(10, ordering=:degrevlex, ground=GF(2^27 + 29)),
            5
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^27+29), cyclic 8",
        result=compute_gb(
            Groebner.cyclicn(8, ordering=:degrevlex, ground=GF(2^27 + 29)),
            5
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), cyclic 8",
        result=compute_gb(Groebner.cyclicn(8, ordering=:degrevlex, ground=GF(2^31 - 1)), 5)
    )
)
push!(
    suite,
    (
        problem_name="groebner, Nemo, GF(2^31-1), cyclic 8",
        result=compute_gb(
            Groebner.cyclicn(8, ordering=:degrevlex, ground=Nemo.GF(2^31 - 1))
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
    minimum(times)
end

# Learn and apply
push!(
    suite,
    (
        problem_name="groebner_apply!, AA, GF(2^31-1), cyclic 7",
        result=learn_and_apply(
            Groebner.cyclicn(7, ordering=:degrevlex, ground=GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner_apply!, Nemo, GF(2^31-1), cyclic 7",
        result=learn_and_apply(
            Groebner.cyclicn(7, ordering=:degrevlex, ground=Nemo.GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner_apply!, AA, GF(2^31-1), katsura 10",
        result=learn_and_apply(
            Groebner.katsuran(10, ordering=:degrevlex, ground=GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner_apply!, AA, GF(2^27+29), katsura 10",
        result=learn_and_apply(
            Groebner.katsuran(10, ordering=:degrevlex, ground=GF(2^27 + 29))
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner_apply!, Nemo, GF(2^31-1), katsura 10",
        result=learn_and_apply(
            Groebner.katsuran(10, ordering=:degrevlex, ground=Nemo.GF(2^31 - 1))
        )
    )
)

# Compute Groebner bases over the rationals
push!(
    suite,
    (
        problem_name="groebner, AA, QQ, katsura 8",
        result=compute_gb(Groebner.katsuran(8, ordering=:degrevlex, ground=QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner, Nemo, QQ, katsura 8",
        result=compute_gb(Groebner.katsuran(8, ordering=:degrevlex, ground=Nemo.QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, QQ, eco 10",
        result=compute_gb(Groebner.eco10(ordering=:degrevlex, ground=QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, QQ, cyclic 7",
        result=compute_gb(Groebner.cyclicn(7, ordering=:degrevlex, ground=QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, QQ, noon 8",
        result=compute_gb(Groebner.noonn(8, ordering=:degrevlex, ground=QQ), 3)
    )
)

function multimodular_gb_problem(nbits; np=AbstractAlgebra)
    R, (x1, x2, x3, x4) =
        polynomial_ring(np.QQ, ["x1", "x2", "x3", "x4"], ordering=:degrevlex)
    nbits_per_prime = 31
    nprimes = max(div(nbits, nbits_per_prime), 1)
    N = prod(map(BigInt, Primes.nextprimes(2^31 - 100, nprimes)))
    @info "Constructing a multi-modular problem" np nbits nbits_per_prime nprimes
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
        problem_name="groebner, AA, QQ, 1000-bit output coeffs",
        result=compute_gb(multimodular_gb_problem(1000, np=AbstractAlgebra))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, QQ, 10000-bit output coeffs",
        result=compute_gb(multimodular_gb_problem(10000, np=AbstractAlgebra))
    )
)
push!(
    suite,
    (
        problem_name="groebner, Nemo, QQ, 1000-bit output coeffs",
        result=compute_gb(multimodular_gb_problem(1000, np=Nemo))
    )
)
push!(
    suite,
    (
        problem_name="groebner, Nemo, QQ, 10000-bit output coeffs",
        result=compute_gb(multimodular_gb_problem(10000, np=Nemo))
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
    minimum(times)
end

# Compute normal forms over integers modulo a prime
push!(
    suite,
    (
        problem_name="normalform, AA, GF(2^31-1), cyclic 7",
        result=compute_normalforms(
            Groebner.cyclicn(7, ordering=:degrevlex, ground=GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="normalform, AA, GF(103), cyclic 8",
        result=compute_normalforms(
            Groebner.cyclicn(8, ordering=:degrevlex, ground=GF(103))
        )
    )
)
push!(
    suite,
    (
        problem_name="normalform, Nemo, GF(103), cyclic 8",
        result=compute_normalforms(
            Groebner.cyclicn(8, ordering=:degrevlex, ground=Nemo.GF(103))
        )
    )
)
push!(
    suite,
    (
        problem_name="normalform, AA, QQ, katsura 9",
        result=compute_normalforms(Groebner.katsuran(9, ordering=:degrevlex, ground=QQ))
    )
)

function dump_results(file, key)
    open(file, "w") do out
        println(out, key)
        for problem in suite
            problem_name = problem.problem_name
            result = problem.result
            println(out, "$problem_name, $result")
        end
    end
end
