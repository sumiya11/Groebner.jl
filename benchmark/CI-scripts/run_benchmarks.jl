using Groebner
using AbstractAlgebra

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

push!(
    suite,
    (
        problem_name="groebner, GF(2^31-1), katsura 5",
        result=compute_gb(Groebner.katsuran(5, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, GF(2^31-1), katsura 6",
        result=compute_gb(Groebner.katsuran(6, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, GF(2^31-1), katsura 8",
        result=compute_gb(Groebner.katsuran(8, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, GF(2^31-1), katsura 10",
        result=compute_gb(Groebner.katsuran(10, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, GF(2^31-1), cyclic 8",
        result=compute_gb(Groebner.cyclicn(8, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, QQ, katsura 8",
        result=compute_gb(Groebner.katsuran(8, ordering=:degrevlex, ground=QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner, QQ, eco 10",
        result=compute_gb(Groebner.eco10(ordering=:degrevlex, ground=QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner, QQ, cyclic 7",
        result=compute_gb(Groebner.cyclicn(7, ordering=:degrevlex, ground=QQ))
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
