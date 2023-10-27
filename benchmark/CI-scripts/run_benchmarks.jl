# This script is not intended to be run directly.
# Instead, see runtests.jl.

using Groebner
using AbstractAlgebra
import Nemo

suite = []

function compute_gb(system)
    times = []
    trials = 7
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
        result=compute_gb(Groebner.katsuran(10, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), cyclic 8",
        result=compute_gb(Groebner.cyclicn(8, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, Nemo, GF(2^31-1), cyclic 8",
        result=compute_gb(
            Groebner.cyclicn(8, ordering=:degrevlex, ground=Nemo.GF(2^31 - 1), np=Nemo)
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
        result=compute_gb(
            Groebner.katsuran(8, ordering=:degrevlex, ground=Nemo.QQ, np=Nemo)
        )
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
            Groebner.cyclicn(8, ordering=:degrevlex, ground=Nemo.GF(103), np=Nemo)
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
