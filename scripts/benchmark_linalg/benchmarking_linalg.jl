using Nemo, Groebner, BenchmarkTools, PrettyTables, Printf

k = Nemo.GF(2^30+3)
systems = [
    ("chandra-9", Groebner.Examples.chandran(9, k=k)),
    ("chandra-10", Groebner.Examples.chandran(10, k=k)),
    ("chandra-11", Groebner.Examples.chandran(11, k=k)),
    ("chandra-12", Groebner.Examples.chandran(12, k=k)),
    ("cyclic-7", Groebner.Examples.cyclicn(7, k=k)),
    ("cyclic-8", Groebner.Examples.cyclicn(8, k=k)),
    ("cyclic-9", Groebner.Examples.cyclicn(9, k=k)),
    ("eco-11", Groebner.Examples.econ(11, k=k)),
    ("eco-12", Groebner.Examples.econ(12, k=k)),
    ("eco-13", Groebner.Examples.econ(13, k=k)),
    ("noon-7", Groebner.Examples.noonn(7, k=k)),
    ("noon-8", Groebner.Examples.noonn(8, k=k)),
    ("noon-9", Groebner.Examples.noonn(9, k=k)),
    ("katsura-9", Groebner.Examples.katsuran(9, k=k)),
    ("katsura-10", Groebner.Examples.katsuran(10, k=k)),
    ("katsura-11", Groebner.Examples.katsuran(11, k=k)),
    ("katsura-12", Groebner.Examples.katsuran(12, k=k)),
    ("hexapod", Groebner.Examples.hexapod(k=k)),
    ("alea6", Groebner.Examples.alea6(k=k)),
    ("hiv2", Groebner.Examples.HIV2(k=k)),
    ("goodwin", Groebner.Examples.Goodwin_with_weights(k=k)),
    ("crn", Groebner.Examples.ChemicalReactionNetwork(k=k)),
    ("jason210", Groebner.Examples.jason210(k=k)),
    ("gametwo2", Groebner.Examples.gametwo2(k=k)),
    ("yang1", Groebner.Examples.yang1(k=k)),
    ("bayes148", Groebner.Examples.bayes148(k=k)),
    ("mayr42", Groebner.Examples.mayr42(k=k)),
]

try
    include("siwr.jl")
    sys = siwr()
    push!(systems, ("siwr", map(f -> map_coefficients(c -> k(numerator(c)) // k(denominator(c)), f), sys)))
catch e @warn "Could not include siwr.jl system: $e"; end

println("Running the following systems: ", map(first, systems))

data = []
for (name, sys) in systems
    @info "Running $name.."
    t1, t2, t3, t4 = 0, 0, 0, 0
    for _ in 1:2
        t1 = @elapsed groebner(sys; threaded=:no)
        t2 = @elapsed trace, _ = groebner_learn(sys)
        t3 = @elapsed groebner_apply!(trace, sys)
        t3 = @elapsed groebner_apply!(trace, sys)
        t4 = @elapsed groebner_apply!(trace, (sys, sys, sys, sys))
        t4 = @elapsed groebner_apply!(trace, (sys, sys, sys, sys))
    end
    push!(data, [name, t1, t2, t3, t4, t1 / t3, 4t1 / t4])
end

matrix = permutedims(reduce(hcat, data))
pretty_table(
    matrix, 
    column_labels=[
        [EmptyCells(2), MultiColumn(3, "Learn & Apply"), MultiColumn(2, "Speedup")], 
        ["Name", "Monte-Carlo", "Learn", "Apply", "Apply 4x", "M.C. / Apply", "4*M.C. / Apply 4x"]
    ],
    title="Timings in seconds",
    formatters=[(v,i,j) -> v isa Number ? @sprintf("%.2f", v) : v],
    table_format = TextTableFormat(borders = text_table_borders__ascii_rounded);  
)

#=
                          Timings in seconds
.------------.-------------.----------------.------------------------.
|            |             | Learn & Apply  |                        |
|       Name | Monte-Carlo |  Learn | Apply | Speedup (M.C. / Apply) |
:------------+-------------+--------+-------+------------------------:
|  chandra-9 |        0.03 |   0.04 |  0.01 |                   2.46 |
| chandra-10 |        0.10 |   0.32 |  0.04 |                   2.52 |
| chandra-11 |        1.18 |   1.75 |  0.28 |                   4.26 |
|   cyclic-7 |        0.14 |   0.12 |  0.03 |                   4.30 |
|   cyclic-8 |        0.86 |   2.83 |  0.86 |                   1.00 |
|     eco-11 |        0.40 |   0.92 |  0.15 |                   2.66 |
|     eco-12 |        1.94 |   6.17 |  1.53 |                   1.27 |
|     noon-7 |        0.13 |   0.31 |  0.05 |                   2.76 |
|     noon-8 |        1.21 |   1.43 |  0.33 |                   3.68 |
|     noon-9 |       11.14 |  11.70 |  3.65 |                   3.05 |
|  katsura-9 |        0.15 |   0.79 |  0.08 |                   1.83 |
| katsura-10 |        0.92 |   3.92 |  0.60 |                   1.53 |
| katsura-11 |        5.31 |  35.84 |  2.72 |                   1.95 |
|    hexapod |        0.00 |   0.01 |  0.00 |                   1.93 |
|      alea6 |        0.09 |   0.20 |  0.21 |                   0.43 |
|       hiv2 |        1.17 |   1.33 |  0.20 |                   5.71 |
|    goodwin |      108.04 | 113.70 | 13.98 |                   7.73 |
|        crn |        0.13 |   0.49 |  0.04 |                   2.86 |
|   jason210 |        2.93 |   3.32 |  0.68 |                   4.34 |
|   gametwo2 |       10.83 |  36.61 |  6.70 |                   1.62 |
|      yang1 |       13.61 |  76.82 |  0.97 |                  14.09 |
|   bayes148 |       29.17 |  58.14 |  0.74 |                  39.52 |
|     mayr42 |       33.47 |  46.58 |  0.18 |                 185.58 |
|       siwr |       17.95 |  16.30 |  0.20 |                  89.91 |
'------------'-------------'--------'-------'------------------------'
=#
