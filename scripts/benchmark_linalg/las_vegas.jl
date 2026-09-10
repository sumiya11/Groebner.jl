using AbstractAlgebra, Groebner, Printf

const LINALG_OPTIONS = (:deterministic, :randomized, :las_vegas)
const SEED = 42
const TASKS = 1

function benchmark_groebner(system, linalg::Symbol, trials::Int)
    times = Vector{Float64}(undef, trials)
    for i in eachindex(times)
        GC.gc()
        times[i] = @elapsed groebner(system; linalg=linalg, seed=SEED, tasks=TASKS)
    end
    minimum(times)
end

function main()
    trials = isempty(ARGS) ? 1 : parse(Int, only(ARGS))

    field = GF(2^30 + 3)
    systems = [
        ("eco-11", Groebner.Examples.econ(11, k=field)),
        ("cholera", Groebner.Examples.Cholera(k=field)),
        ("goodwin", Groebner.Examples.Goodwin_with_weights(k=field)),
        ("cyclic-9", Groebner.Examples.cyclicn(9, k=field)),
        ("noon-9", Groebner.Examples.noonn(9, k=field)),
        ("katsura-11", Groebner.Examples.katsuran(11, k=field)),
    ]

    println("Groebner linear-algebra benchmark")
    println("field: GF(2^30 + 3), trials: $trials, tasks: $TASKS, seed: $SEED")
    println("Each entry is the best runtime over all trials, after warmup.")
    println()
    @printf(
        "%-12s %15s %15s %15s %15s\n",
        "system",
        "deterministic",
        "randomized",
        "las-vegas",
        "LV/det."
    )

    for (name, system) in systems
        bases = map(
            linalg -> groebner(system; linalg=linalg, seed=SEED, tasks=TASKS),
            LINALG_OPTIONS
        )
        all(basis -> basis == first(bases), bases) ||
            error("linear algebra backends disagree on $name")

        deterministic = benchmark_groebner(system, :deterministic, trials)
        randomized = benchmark_groebner(system, :randomized, trials)
        las_vegas = benchmark_groebner(system, :las_vegas, trials)
        @printf(
            "%-12s %12.3f s %12.3f s %12.3f s %15.2f\n",
            name,
            deterministic,
            randomized,
            las_vegas,
            las_vegas / deterministic
        )
    end
end

main()
