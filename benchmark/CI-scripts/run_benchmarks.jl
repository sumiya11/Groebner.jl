using Groebner
using AbstractAlgebra

suite = []

function compute_gb(system)
    times = []
    trials = 5
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
        problem_name="finite field, katsura 5",
        result=compute_gb(Groebner.katsuran(5, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="finite field, katsura 6",
        result=compute_gb(Groebner.katsuran(6, ordering=:degrevlex, ground=GF(2^31 - 1)))
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
