using AbstractAlgebra, Groebner, BenchmarkTools, PrettyTables, Printf, TimerOutputs

k = AbstractAlgebra.QQ
system_solving = [
    ("chandra-11", Groebner.Examples.chandran(11, k=k)),
    ("chandra-12", Groebner.Examples.chandran(12, k=k)),
    ("cyclic-8", Groebner.Examples.cyclicn(8, k=k)),
    # ("cyclic-9", Groebner.Examples.cyclicn(9, k=k)),
    ("katsura-10", Groebner.Examples.katsuran(10, k=k)),
    ("katsura-11", Groebner.Examples.katsuran(11, k=k)),
    ("reimer-8", Groebner.Examples.reimern(8, k=k)),
    ("ipp", Groebner.Examples.ipp(k=k)),
]
sian = [
    # ("cholera", Groebner.Examples.Cholera(k=k)),
    ("hiv2", Groebner.Examples.HIV2(k=k)),
    # ("goodwin", Groebner.Examples.Goodwin_with_weights(k=k)),
    # ("crn", Groebner.Examples.ChemicalReactionNetwork(k=k)),
]
other = [
    ("alea6", Groebner.Examples.alea6(k=k)),
    # ("jason210", Groebner.Examples.jason210(k=k)),
    # ("gametwo2", Groebner.Examples.gametwo2(k=k)),
    # ("yang1", Groebner.Examples.yang1(k=k)),
    # ("bayes148", Groebner.Examples.bayes148(k=k)),
    # ("mayr42", Groebner.Examples.mayr42(k=k))
]
# random = sort(reduce(vcat, [
#     ("rand-$(n)-$d", randsys(n, d))
#     for n in 4:6, d in 4:6
# ]), by=first)
# random = sort(reduce(vcat, [
#     ("rand-$(n)-$d", randsys(n, d))
#     for n in 8:14, d in 2:2
# ]), by=first)
random = []
systems = vcat(system_solving, sian, other, random)

println("Running the following systems: ", map(first, systems))

TimerOutputs.disable_timer!(Groebner._TIMER)
timers = []
for (name, sys) in systems
    @info "Running $name.."

    for kws in [
        (tasks=1,),
        (tasks=2,),
        (tasks=4,),
        (tasks=8,),
    ]
        @info "Options:" kws
        for t in 1:2
            TimerOutputs.enable_timer!(Groebner._TIMER2); reset_timer!(Groebner._TIMER2);
            groebner(sys; kws...);
            if t == 2 show(Groebner._TIMER2, allocations=false); println() end
        end
        push!(timers, [name, kws, copy(Groebner._TIMER2)])
    end    
end

labels = ["Name", "Options", "Guess Prime", "Learn", "Apply", "CRT", "RatRec", "Check", "Total (recorded)"]
data = []
for entry in timers
    name = entry[1]
    kws = entry[2]
    timer = entry[3]
    time_guess_prime = TimerOutputs.time(timer["_groebner_learn_and_apply"]["_groebner_guess_lucky_prime"])
    time_learn = TimerOutputs.time(timer["_groebner_learn_and_apply"]["f4_learn!"])
    time_apply = TimerOutputs.time(timer["_groebner_learn_and_apply"]["f4_apply!"])
    time_crt = TimerOutputs.time(timer["_groebner_learn_and_apply"]["crt_vec_full!"]) + 
               TimerOutputs.time(timer["_groebner_learn_and_apply"]["crt_vec_partial!"])
    time_ratrec = TimerOutputs.time(timer["_groebner_learn_and_apply"]["ratrec_vec_full!"]) + 
                  TimerOutputs.time(timer["_groebner_learn_and_apply"]["ratrec_vec_partial!"])
    # time_io = TimerOutputs.time(timer["io_convert_polynomials_to_ir"]) +
    #           TimerOutputs.time(timer["io_convert_ir_to_polynomials"]) +
    #           TimerOutputs.time(timer["ir_convert_internal_to_ir"]) +
    #             TimerOutputs.time(timer["ir_convert_ir_to_internal"])
    time_check = TimerOutputs.time(timer["_groebner_learn_and_apply"]["modular_lift_check!"])
    time_total_recorded = TimerOutputs.tottime(timer)
    push!(data, [name, kws, time_guess_prime, time_learn, time_apply, time_crt, time_ratrec, time_check, time_total_recorded])
end

matrix = permutedims(reduce(hcat, data))

threshold = 0.2
hl5 = TextHighlighter((v,i,j) -> (total = v[i,end]; j != size(matrix, 2) && v[i,j] isa Number && threshold < v[i,j] / total < 2*threshold), crayon"bg:(255,220,220)")
hl6 = TextHighlighter((v,i,j) -> (total = v[i,end]; j != size(matrix, 2) && v[i,j] isa Number && 2*threshold < v[i,j] / total), crayon"bg:(255,170,170)")

pretty_table(
    matrix, 
    column_labels=labels,
    title="Timings in seconds",
    formatters=[
        # (v,i,j) -> v isa Number ? (j != size(matrix, 2) ? "$(round(Int, 100*v / matrix[i,end]))%" : @sprintf("%.2f", v / 1e9)) : v,
        (v,i,j) -> v isa Number ? (j != size(matrix, 2) ? @sprintf("%.2f", v / 1e9) : @sprintf("%.2f", v / 1e9)) : v,
    ],
    row_group_labels = [1+4*(i-1) => systems[i][1] for i in 1:length(systems)],
    highlighters  = [hl5, hl6],
    table_format = TextTableFormat(borders = text_table_borders__ascii_rounded),
    # summary_row_labels = ["Total"],
    # summary_rows = [(data, i) -> i > 1 ? @sprintf("%d%%", round(Int, 100*sum(data[:, i]) / sum(data[:, end]))) : ""],
    fit_table_in_display_vertically = false,
)

#=
                                          Timings in seconds
.------------.--------------.-------------.-------.-------.-------.--------.-------.------------------.
|       Name |      Options | Guess Prime | Learn | Apply |   CRT | RatRec | Check | Total (recorded) |
:------------'--------------'-------------'-------'-------'-------'--------'-------'------------------:
| chandra-11                                                                                          |
:------------.--------------.-------------.-------.-------.-------.--------.-------.------------------:
| chandra-11 | (tasks = 1,) |        0.83 |  0.97 |  1.59 |  2.56 |   1.77 |  0.50 |             8.42 |
| chandra-11 | (tasks = 2,) |        0.81 |  0.97 |  1.25 |  1.69 |   1.29 |  0.54 |             6.80 |
| chandra-11 | (tasks = 4,) |        0.93 |  1.06 |  0.93 |  1.04 |   0.78 |  0.58 |             5.49 |
| chandra-11 | (tasks = 8,) |        0.71 |  1.01 |  0.64 |  0.81 |   0.81 |  0.50 |             4.85 |
:------------'--------------'-------------'-------'-------'-------'--------'-------'------------------:
| chandra-12                                                                                          |
:------------.--------------.-------------.-------.-------.-------.--------.-------.------------------:
| chandra-12 | (tasks = 1,) |        3.63 |  5.69 |  8.75 | 13.65 |   7.59 |  1.84 |            41.56 |
| chandra-12 | (tasks = 2,) |        3.73 |  5.80 |  5.75 |  8.03 |   4.37 |  1.86 |            29.90 |
| chandra-12 | (tasks = 4,) |        3.67 |  5.76 |  3.60 |  4.85 |   3.06 |  1.69 |            23.10 |
| chandra-12 | (tasks = 8,) |        3.82 |  5.74 |  2.69 |  3.10 |   3.43 |  2.17 |            21.50 |
:------------'--------------'-------------'-------'-------'-------'--------'-------'------------------:
| cyclic-8                                                                                            |
:------------.--------------.-------------.-------.-------.-------.--------.-------.------------------:
|   cyclic-8 | (tasks = 1,) |        1.59 |  2.59 |  5.72 |  0.44 |   0.84 |  0.17 |            11.41 |
|   cyclic-8 | (tasks = 2,) |        1.55 |  2.56 |  3.21 |  0.18 |   0.39 |  0.13 |             8.07 |
|   cyclic-8 | (tasks = 4,) |        1.51 |  2.55 |  1.75 |  0.19 |   0.25 |  0.10 |             6.39 |
|   cyclic-8 | (tasks = 8,) |        1.59 |  2.63 |  3.71 |  0.34 |   0.35 |  0.19 |             8.91 |
:------------'--------------'-------------'-------'-------'-------'--------'-------'------------------:
| katsura-10                                                                                          |
:------------.--------------.-------------.-------.-------.-------.--------.-------.------------------:
| katsura-10 | (tasks = 1,) |        1.21 |  3.65 |  3.70 |  1.38 |   1.96 |  0.63 |            12.70 |
| katsura-10 | (tasks = 2,) |        1.24 |  3.62 |  1.89 |  0.71 |   1.21 |  0.61 |             9.47 |
| katsura-10 | (tasks = 4,) |        1.13 |  3.59 |  1.71 |  0.89 |   0.89 |  0.59 |             9.04 |
| katsura-10 | (tasks = 8,) |        1.13 |  3.62 |  1.92 |  0.86 |   0.96 |  0.74 |             9.51 |
:------------'--------------'-------------'-------'-------'-------'--------'-------'------------------:
| katsura-11                                                                                          |
:------------.--------------.-------------.-------.-------.-------.--------.-------.------------------:
| katsura-11 | (tasks = 1,) |        8.82 | 40.77 | 39.15 |  7.92 |   9.81 |  1.75 |           109.26 |
| katsura-11 | (tasks = 2,) |        8.63 | 40.68 | 26.72 |  4.71 |   5.61 |  1.92 |            88.94 |
| katsura-11 | (tasks = 4,) |        8.60 | 40.53 | 16.80 |  4.28 |   4.51 |  2.07 |            77.56 |
| katsura-11 | (tasks = 8,) |        8.87 | 41.56 | 13.61 |  3.47 |   3.85 |  1.97 |            74.36 |
:------------'--------------'-------------'-------'-------'-------'--------'-------'------------------:
| reimer-8                                                                                            |
:------------.--------------.-------------.-------.-------.-------.--------.-------.------------------:
|   reimer-8 | (tasks = 1,) |       22.53 | 37.70 | 33.43 | 10.97 |   9.01 |  3.92 |           118.81 |
|   reimer-8 | (tasks = 2,) |       22.70 | 38.24 | 18.98 |  6.32 |   4.86 |  4.19 |            96.02 |
|   reimer-8 | (tasks = 4,) |       23.04 | 39.06 | 12.04 |  5.48 |   3.94 |  3.75 |            88.30 |
|   reimer-8 | (tasks = 8,) |       23.97 | 40.42 | 12.75 |  4.85 |   3.57 |  4.09 |            91.02 |
:------------'--------------'-------------'-------'-------'-------'--------'-------'------------------:
| ipp                                                                                                 |
:------------.--------------.-------------.-------.-------.-------.--------.-------.------------------:
|        ipp | (tasks = 1,) |        0.01 |  0.01 |  1.30 | 22.44 |   5.16 |  0.04 |            28.99 |
|        ipp | (tasks = 2,) |        0.01 |  0.01 |  0.71 | 10.38 |   3.23 |  0.05 |            14.42 |
|        ipp | (tasks = 4,) |        0.01 |  0.01 |  0.40 |  8.39 |   2.53 |  0.04 |            11.41 |
|        ipp | (tasks = 8,) |        0.01 |  0.01 |  0.21 |  6.59 |   2.05 |  0.04 |             8.94 |
:------------'--------------'-------------'-------'-------'-------'--------'-------'------------------:
| hiv2                                                                                                |
:------------.--------------.-------------.-------.-------.-------.--------.-------.------------------:
|       hiv2 | (tasks = 1,) |        2.34 |  1.35 |  8.48 |  0.01 |   0.01 |  0.00 |            12.25 |
|       hiv2 | (tasks = 2,) |        2.33 |  1.35 |  4.47 |  0.01 |   0.00 |  0.00 |             8.25 |
|       hiv2 | (tasks = 4,) |        2.36 |  1.37 |  2.54 |  0.01 |   0.01 |  0.00 |             6.47 |
|       hiv2 | (tasks = 8,) |        2.30 |  1.48 |  6.21 |  0.03 |   0.00 |  0.03 |            10.40 |
:------------'--------------'-------------'-------'-------'-------'--------'-------'------------------:
| alea6                                                                                               |
:------------.--------------.-------------.-------.-------.-------.--------.-------.------------------:
|      alea6 | (tasks = 1,) |        0.16 |  0.23 |  5.39 |  7.26 |   6.19 |  0.29 |            19.56 |
|      alea6 | (tasks = 2,) |        0.21 |  0.18 |  3.10 |  4.75 |   3.61 |  0.25 |            12.16 |
|      alea6 | (tasks = 4,) |        0.15 |  0.18 |  1.91 |  2.38 |   2.14 |  0.29 |             7.16 |
|      alea6 | (tasks = 8,) |        0.21 |  0.17 |  1.53 |  1.79 |   1.44 |  0.27 |             5.49 |
'------------'--------------'-------------'-------'-------'-------'--------'-------'------------------'

=#
