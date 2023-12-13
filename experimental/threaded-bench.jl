using BenchmarkTools, AbstractAlgebra, PrettyTables, Groebner
using Base.Threads

@info "Using $(nthreads()) Julia threads"

Groebner.logging_enabled() = false
Groebner.invariants_enabled() = false

p = 2^30 + 3
systems = [
    ("kat5", Groebner.katsuran(5, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("kat6", Groebner.katsuran(6, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("kat7", Groebner.katsuran(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("kat8", Groebner.katsuran(8, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("kat9", Groebner.katsuran(9, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("kat10", Groebner.katsuran(10, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("hen5", Groebner.henrion5(ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("reim4", Groebner.reimern(4, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("reim5", Groebner.reimern(5, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("cyc4", Groebner.cyclicn(4, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("cyc5", Groebner.cyclicn(5, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("cyc6", Groebner.cyclicn(6, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("cyc7", Groebner.cyclicn(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("cyc8", Groebner.cyclicn(8, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("eco10", Groebner.eco10(ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("eco11", Groebner.eco11(ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("eco12", Groebner.eco12(ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("noon4", Groebner.noonn(4, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("noon5", Groebner.noonn(5, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("noon6", Groebner.noonn(6, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("noon7", Groebner.noonn(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("noon8", Groebner.noonn(8, ordering=:degrevlex, ground=AbstractAlgebra.GF(p)))
]

table = Matrix{Any}(undef, (length(systems), 3))

function benchmark_system(system, trials=5; kwargs...)
    timings = []
    GC.gc()
    gb = nothing
    for _ in 1:trials
        time = @elapsed gb = Groebner.groebner(system; kwargs...)
        push!(timings, time)
    end
    gb, minimum(timings)
end

for (i, (name, s)) in enumerate(systems)
    @info """
    $name:
    lock_free / compare_and_swap / single_thread"""
    gb1, t1 =
        benchmark_system(s; linalg=:deterministic, threaded=:yes, threadalgo=:lock_free)
    gb2, t2 = benchmark_system(
        s;
        linalg=:deterministic,
        threaded=:yes,
        threadalgo=:compare_and_swap
    )
    gb3, t3 = benchmark_system(s; linalg=:deterministic, threaded=:no)

    (t1, t2, t3) = map(t -> BenchmarkTools.prettytime(t * 1e9), (t1, t2, t3))
    table[i, :] .= (t1, t2, t3)
    println("$t1 / $t2 / $t3")

    @assert gb1 == gb2 == gb3
end

pretty_table(
    table,
    header=["lock_free", "compare_and_swap", "single_thread"],
    tf=tf_markdown,
    row_labels=map(first, systems)
)
