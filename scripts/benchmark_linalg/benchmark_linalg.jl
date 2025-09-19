using AbstractAlgebra, Groebner, BenchmarkTools, PrettyTables, Printf, TimerOutputs

include("gizmos.jl")

k = AbstractAlgebra.GF(2^30+3)
system_solving = [
    ("chandra-9", Groebner.Examples.chandran(9, k=k)),
    ("chandra-10", Groebner.Examples.chandran(10, k=k)),
    ("chandra-11", Groebner.Examples.chandran(11, k=k)),
    ("chandra-12", Groebner.Examples.chandran(12, k=k)),
    ("chandra-13", Groebner.Examples.chandran(13, k=k)),
    ("cyclic-7", Groebner.Examples.cyclicn(7, k=k)),
    ("cyclic-8", Groebner.Examples.cyclicn(8, k=k)),
    # ("cyclic-9", Groebner.Examples.cyclicn(9, k=k)),
    ("eco-11", Groebner.Examples.econ(11, k=k)),
    ("eco-12", Groebner.Examples.econ(12, k=k)),
    ("eco-13", Groebner.Examples.econ(13, k=k)),
#    ("eco-14", Groebner.Examples.econ(14, k=k)),
    ("noon-7", Groebner.Examples.noonn(7, k=k)),
    ("noon-8", Groebner.Examples.noonn(8, k=k)),
    ("noon-9", Groebner.Examples.noonn(9, k=k)),
    ("noon-10", Groebner.Examples.noonn(10, k=k)),
    ("henrion-6", Groebner.Examples.henrion6(k=k)),
    ("henrion-7", Groebner.Examples.henrion7(k=k)),
#    ("henrion-8", Groebner.Examples.henrion8(k=k)),
    ("katsura-9", Groebner.Examples.katsuran(9, k=k)),
    ("katsura-10", Groebner.Examples.katsuran(10, k=k)),
    ("katsura-11", Groebner.Examples.katsuran(11, k=k)),
    ("katsura-12", Groebner.Examples.katsuran(12, k=k)),
    ("reimer-6", Groebner.Examples.reimern(6, k=k)),
    ("reimer-7", Groebner.Examples.reimern(7, k=k)),
    ("reimer-8", Groebner.Examples.reimern(8, k=k)),
    ("hexapod", Groebner.Examples.hexapod(k=k)),
]
sian = [
    ("cholera", Groebner.Examples.Cholera(k=k)),
    ("hiv2", Groebner.Examples.HIV2(k=k)),
    ("goodwin", Groebner.Examples.Goodwin_with_weights(k=k)),
    ("crn", Groebner.Examples.ChemicalReactionNetwork(k=k)),

]
other = [
    ("alea6", Groebner.Examples.alea6(k=k)),
    ("jason210", Groebner.Examples.jason210(k=k)),
    ("gametwo2", Groebner.Examples.gametwo2(k=k)),
    ("yang1", Groebner.Examples.yang1(k=k)),
    ("bayes148", Groebner.Examples.bayes148(k=k)),
    ("mayr42", Groebner.Examples.mayr42(k=k))
]
random = sort(reduce(vcat, [
    ("rand-$(n)-$d", randsys(n, d))
    for n in 4:5, d in 4:5
]), by=first)
# random = sort(reduce(vcat, [
#     ("rand-$(n)-$d", randsys(n, d))
#     for n in 8:14, d in 2:2
# ]), by=first)
systems = vcat(system_solving, sian, other, random)
# systems = vcat(random)

println("Running the following systems: ", map(first, systems))

timers = []
for (name, sys) in systems
    @info "Running $name.."

    for kws in [
        1,
        2,
        3
    ]
        opt = nothing
        @info "Options:" kws
        for t in 1:2
            Groebner._LINALG_MOD_P[] = Groebner._LINALG_ADDMUL_MOD_P[] = Groebner._LINALG_MUL_MOD_P[] = Groebner._LINALG_REDUCE_ROW[] = (0, 0)
            Groebner._LINALG_REDUCER_ROWS[] = 0
            TimerOutputs.enable_timer!(Groebner._TIMER); reset_timer!(Groebner._TIMER);
            if kws == 2 
              groebner(sys; linalg=:randomized);
              opt = (linalg=:randomized,)
      elseif kws == 1
              groebner(sys; linalg=:deterministic);
              opt = (linalg=:deterministic,)
      elseif kws == 3
              trace, _ = groebner_learn(sys);
              Groebner._LINALG_MOD_P[] = Groebner._LINALG_ADDMUL_MOD_P[] = Groebner._LINALG_MUL_MOD_P[] = Groebner._LINALG_REDUCE_ROW[] = (0, 0)
              Groebner._LINALG_REDUCER_ROWS[] = 0
              TimerOutputs.enable_timer!(Groebner._TIMER); reset_timer!(Groebner._TIMER);
              groebner_apply!(trace, sys);
              opt = (:groebner_apply,)
      else
        error("beda")
      end
            if t == 2 show(Groebner._TIMER, allocations=false); println() end
        end
        push!(timers, [name, opt, copy(Groebner._TIMER), Groebner._LINALG_MOD_P[], Groebner._LINALG_ADDMUL_MOD_P[], Groebner._LINALG_MUL_MOD_P[], Groebner._LINALG_REDUCE_ROW[], Groebner._LINALG_REDUCER_ROWS[]])
    end    
end

#=
Base.getindex(::Nothing,key::String) = nothing
Base.getindex(timer::TimerOutput,key::String) = haskey(timer, key) ? timer.inner_timers[key] : nothing
TimerOutputs.time(::Nothing) = 0
=#

labels = ["Name", "Options", "Select Pairs", "Symb Prepr", "Reduction", "Update", "Autoreduce", "IO", "Total"]
labels2 = ["Name", "Options", "% p", "addmul", "mul", "rows reduced", "rows reducers"]
data = []
data2 = []
for entry in timers
    name = entry[1]
    kws = entry[2]
    timer = entry[3]
    stats = entry[4:end]
    if kws == (:groebner_apply,)
      time_select_pairs = 0
      time_symb_preprc = TimerOutputs.time(timer["f4_apply!"]["f4_symbolic_preprocessing!"])
      time_reduction = TimerOutputs.time(timer["f4_apply!"]["f4_reduction_apply!"])
      time_update = TimerOutputs.time(timer["f4_apply!"]["basis_update!"])
      time_autoreduce = TimerOutputs.time(timer["f4_apply!"]["f4_autoreduce_apply!"])
      time_io = TimerOutputs.time(timer["io_convert_polynomials_to_ir"]) +
              TimerOutputs.time(timer["io_convert_ir_to_polynomials"])
    else
      time_select_pairs = TimerOutputs.time(timer["f4!"]["f4_select_critical_pairs!"])
      time_symb_preprc = TimerOutputs.time(timer["f4!"]["f4_symbolic_preprocessing!"])
      time_reduction = TimerOutputs.time(timer["f4!"]["f4_reduction!"])
      time_update = TimerOutputs.time(timer["f4!"]["f4_update!"])
      time_autoreduce = TimerOutputs.time(timer["f4!"]["f4_autoreduce!"])
      time_io = TimerOutputs.time(timer["io_convert_polynomials_to_ir"]) +
              TimerOutputs.time(timer["io_convert_ir_to_polynomials"]) +
              TimerOutputs.time(timer["ir_convert_internal_to_ir"]) +
                TimerOutputs.time(timer["ir_convert_ir_to_internal"])
    end
    time_total = TimerOutputs.tottime(timer)
    push!(data, [name, kws, time_select_pairs, time_symb_preprc, time_reduction, time_update, time_autoreduce, time_io, time_total])
    push!(data2, [name, kws, stats...])
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
        (v,i,j) -> v isa Number ? @sprintf("%.2f", v / 1e9) : v, 
    ],
    # row_group_labels = [1 => "System solving", length(system_solving)+1 => "SIAN", length(system_solving)+length(sian)+1 => "Other"],
    highlighters  = [hl5, hl6],
    table_format = TextTableFormat(borders = text_table_borders__ascii_rounded),
    summary_row_labels = ["Total"],
    summary_rows = [(data, i) -> i > 2 ? @sprintf("%.2f", sum(data[:, i]) / 1e9) : ""],
    fit_table_in_display_vertically = false,
)

matrix2 = permutedims(reduce(hcat, data2))

pretty_table(
    matrix2, 
    column_labels=labels2,
    title="Count linalg ops",
    # formatters=[
    #     (v,i,j) -> v isa Number ? @sprintf("%.2f", v / 1e9) : v, 
    # ],
    # row_group_labels = [1 => "System solving", length(system_solving)+1 => "SIAN", length(system_solving)+length(sian)+1 => "Other"],
    # highlighters  = [hl_select, hl_symb_prepr, hl_update, hl_io],
    table_format = TextTableFormat(borders = text_table_borders__ascii_rounded),
    summary_row_labels = ["Total"],
    summary_rows = [(data, i) -> i > 2 ? reduce(.+, data[:, i]) : ""],
    fit_table_in_display_vertically = false,
)

#=
                                                        Timings in seconds
.-------.------------.----------------------------.--------------.------------.-----------.--------.------------.-------.---------.
|       |       Name |                    Options | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |    IO |   Total |
:-------+------------+----------------------------+--------------+------------+-----------+--------+------------+-------+---------:
|       |  chandra-9 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |    0.02 |
|       |  chandra-9 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.01 |
|       |  chandra-9 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.00 |
|       | chandra-10 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.07 |   0.01 |       0.00 |  0.00 |    0.09 |
|       | chandra-10 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.02 |   0.01 |       0.00 |  0.00 |    0.04 |
|       | chandra-10 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.01 |
|       | chandra-11 | (linalg = :deterministic,) |         0.01 |       0.05 |      0.44 |   0.03 |       0.01 |  0.10 |    0.64 |
|       | chandra-11 |    (linalg = :randomized,) |         0.02 |       0.02 |      0.12 |   0.03 |       0.01 |  0.02 |    0.21 |
|       | chandra-11 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.02 |
|       | chandra-12 | (linalg = :deterministic,) |         0.05 |       0.13 |      2.76 |   0.19 |       0.03 |  0.10 |    3.25 |
|       | chandra-12 |    (linalg = :randomized,) |         0.05 |       0.22 |      0.67 |   0.10 |       0.02 |  0.08 |    1.15 |
|       | chandra-12 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.05 |    0.11 |
|       | chandra-13 | (linalg = :deterministic,) |         0.19 |       0.63 |     18.24 |   0.43 |       0.18 |  0.62 |   20.29 |
|       | chandra-13 |    (linalg = :randomized,) |         0.26 |       0.49 |      4.19 |   0.58 |       0.11 |  0.59 |    6.22 |
|       | chandra-13 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.28 |    0.64 |
|       |   cyclic-7 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.06 |   0.01 |       0.00 |  0.00 |    0.07 |
|       |   cyclic-7 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.02 |   0.01 |       0.00 |  0.00 |    0.03 |
|       |   cyclic-7 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.01 |
|       |   cyclic-8 | (linalg = :deterministic,) |         0.03 |       0.03 |      1.61 |   0.03 |       0.00 |  0.00 |    1.71 |
|       |   cyclic-8 |    (linalg = :randomized,) |         0.03 |       0.03 |      0.38 |   0.04 |       0.00 |  0.00 |    0.48 |
|       |   cyclic-8 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.25 |
|       |     eco-11 | (linalg = :deterministic,) |         0.01 |       0.02 |      0.43 |   0.01 |       0.02 |  0.00 |    0.49 |
|       |     eco-11 |    (linalg = :randomized,) |         0.01 |       0.02 |      0.07 |   0.01 |       0.02 |  0.01 |    0.15 |
|       |     eco-11 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.01 |    0.07 |
|       |     eco-12 | (linalg = :deterministic,) |         0.03 |       0.09 |      3.33 |   0.11 |       0.14 |  0.05 |    3.74 |
|       |     eco-12 |    (linalg = :randomized,) |         0.03 |       0.18 |      0.43 |   0.04 |       0.14 |  0.03 |    0.85 |
|       |     eco-12 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.43 |
|       |     eco-13 | (linalg = :deterministic,) |         0.21 |       0.48 |     30.02 |   0.10 |       0.04 |  0.31 |   31.15 |
|       |     eco-13 |    (linalg = :randomized,) |         0.10 |       0.60 |      2.71 |   0.08 |       0.06 |  0.41 |    3.96 |
|       |     eco-13 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.18 |    1.69 |
|       |     eco-14 | (linalg = :deterministic,) |         0.72 |       2.46 |    261.37 |   0.36 |       9.92 |  1.22 |  276.06 |
|       |     eco-14 |    (linalg = :randomized,) |         0.53 |       2.52 |     21.09 |   0.30 |       9.93 |  0.85 |   35.21 |
|       |     eco-14 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.65 |   21.24 |
|       |     noon-7 | (linalg = :deterministic,) |         0.01 |       0.02 |      0.03 |   0.01 |       0.00 |  0.00 |    0.07 |
|       |     noon-7 |    (linalg = :randomized,) |         0.00 |       0.02 |      0.03 |   0.01 |       0.00 |  0.00 |    0.06 |
|       |     noon-7 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.01 |
|       |     noon-8 | (linalg = :deterministic,) |         0.03 |       0.13 |      0.46 |   0.05 |       0.02 |  0.04 |    0.73 |
|       |     noon-8 |    (linalg = :randomized,) |         0.05 |       0.19 |      0.27 |   0.05 |       0.02 |  0.03 |    0.61 |
|       |     noon-8 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.03 |    0.11 |
|       |     noon-9 | (linalg = :deterministic,) |         0.17 |       0.63 |      2.92 |   0.38 |       0.14 |  0.53 |    4.77 |
|       |     noon-9 |    (linalg = :randomized,) |         0.14 |       0.64 |      2.88 |   0.38 |       0.14 |  0.42 |    4.60 |
|       |     noon-9 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.02 |    0.67 |
|       |    noon-10 | (linalg = :deterministic,) |         1.13 |       5.54 |     27.96 |   3.42 |       1.19 |  2.48 |   41.72 |
|       |    noon-10 |    (linalg = :randomized,) |         1.05 |       5.94 |     31.25 |   3.42 |       1.00 |  2.11 |   44.77 |
|       |    noon-10 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  1.13 |    5.07 |
|       |  henrion-6 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |    0.02 |
|       |  henrion-6 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |    0.01 |
|       |  henrion-6 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.01 |
|       |  henrion-7 | (linalg = :deterministic,) |         0.03 |       0.28 |      1.89 |   0.01 |       0.01 |  0.09 |    2.31 |
|       |  henrion-7 |    (linalg = :randomized,) |         0.09 |       0.28 |      0.52 |   0.01 |       0.01 |  0.10 |    1.01 |
|       |  henrion-7 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.01 |    0.39 |
|       |  henrion-8 | (linalg = :deterministic,) |         2.43 |      20.51 |    649.74 |   0.19 |       0.88 |  7.62 |  681.38 |
|       |  henrion-8 |    (linalg = :randomized,) |         1.95 |      19.21 |    118.73 |   0.19 |       0.62 |  6.61 |  147.30 |
|       |  henrion-8 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.94 |   78.77 |
|       |  katsura-9 | (linalg = :deterministic,) |         0.00 |       0.01 |      0.26 |   0.00 |       0.00 |  0.00 |    0.28 |
|       |  katsura-9 |    (linalg = :randomized,) |         0.00 |       0.01 |      0.04 |   0.00 |       0.00 |  0.00 |    0.06 |
|       |  katsura-9 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.11 |
|       | katsura-10 | (linalg = :deterministic,) |         0.02 |       0.03 |      2.15 |   0.01 |       0.10 |  0.03 |    2.33 |
|       | katsura-10 |    (linalg = :randomized,) |         0.02 |       0.06 |      0.25 |   0.01 |       0.01 |  0.11 |    0.45 |
|       | katsura-10 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.02 |    0.15 |
|       | katsura-11 | (linalg = :deterministic,) |         0.07 |       0.14 |     17.78 |   0.04 |       0.02 |  0.33 |   18.38 |
|       | katsura-11 |    (linalg = :randomized,) |         0.09 |       0.21 |      1.79 |   0.04 |       0.02 |  0.16 |    2.31 |
|       | katsura-11 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.04 |    0.90 |
|       | katsura-12 | (linalg = :deterministic,) |         0.58 |       0.77 |    160.65 |   0.16 |       0.08 |  1.08 |  163.32 |
|       | katsura-12 |    (linalg = :randomized,) |         0.50 |       0.55 |     12.58 |   0.15 |       0.08 |  1.20 |   15.07 |
|       | katsura-12 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.45 |    7.29 |
|       |   reimer-6 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.02 |   0.00 |       0.00 |  0.00 |    0.03 |
|       |   reimer-6 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |    0.02 |
|       |   reimer-6 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.00 |
|       |   reimer-7 | (linalg = :deterministic,) |         0.01 |       0.09 |      0.47 |   0.00 |       0.00 |  0.01 |    0.58 |
|       |   reimer-7 |    (linalg = :randomized,) |         0.01 |       0.10 |      0.12 |   0.00 |       0.00 |  0.01 |    0.24 |
|       |   reimer-7 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.06 |
|       |   reimer-8 | (linalg = :deterministic,) |         0.24 |       2.20 |     16.59 |   0.02 |       0.03 |  0.15 |   19.23 |
|       |   reimer-8 |    (linalg = :randomized,) |         0.11 |       2.08 |      3.08 |   0.02 |       0.03 |  0.34 |    5.67 |
|       |   reimer-8 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.02 |    0.98 |
|       |    hexapod | (linalg = :deterministic,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.00 |
|       |    hexapod |    (linalg = :randomized,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.00 |
|       |    hexapod |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.00 |
|       |    cholera | (linalg = :deterministic,) |         0.99 |       0.45 |     42.17 |   1.51 |       0.00 |  0.00 |   45.12 |
|       |    cholera |    (linalg = :randomized,) |         0.96 |       0.33 |     13.39 |   1.57 |       0.00 |  0.00 |   16.25 |
|       |    cholera |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    8.74 |
|       |       hiv2 | (linalg = :deterministic,) |         0.05 |       0.14 |      0.44 |   0.33 |       0.00 |  0.00 |    0.96 |
|       |       hiv2 |    (linalg = :randomized,) |         0.05 |       0.17 |      0.32 |   0.39 |       0.00 |  0.00 |    0.94 |
|       |       hiv2 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.12 |
|       |    goodwin | (linalg = :deterministic,) |         4.06 |       1.94 |     44.79 |  14.72 |       0.00 |  0.00 |   65.51 |
|       |    goodwin |    (linalg = :randomized,) |         3.96 |       2.24 |     41.67 |  14.34 |       0.00 |  0.00 |   62.21 |
|       |    goodwin |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    8.64 |
|       |        crn | (linalg = :deterministic,) |         0.01 |       0.02 |      0.28 |   0.01 |       0.00 |  0.00 |    0.32 |
|       |        crn |    (linalg = :randomized,) |         0.01 |       0.02 |      0.05 |   0.01 |       0.00 |  0.00 |    0.09 |
|       |        crn |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.03 |
|       |      alea6 | (linalg = :deterministic,) |         0.00 |       0.06 |      0.09 |   0.00 |       0.00 |  0.00 |    0.16 |
|       |      alea6 |    (linalg = :randomized,) |         0.00 |       0.01 |      0.15 |   0.00 |       0.00 |  0.01 |    0.17 |
|       |      alea6 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.02 |
|       |   jason210 | (linalg = :deterministic,) |         0.08 |       0.75 |      0.66 |   0.04 |       0.17 |  0.04 |    1.74 |
|       |   jason210 |    (linalg = :randomized,) |         0.08 |       0.67 |      0.76 |   0.04 |       0.26 |  0.01 |    1.82 |
|       |   jason210 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.02 |    0.23 |
|       |   gametwo2 | (linalg = :deterministic,) |         0.10 |       0.19 |     20.71 |   0.11 |       0.38 |  0.47 |   21.95 |
|       |   gametwo2 |    (linalg = :randomized,) |         0.23 |       0.14 |      4.21 |   0.10 |       0.40 |  0.53 |    5.61 |
|       |   gametwo2 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.28 |    3.62 |
|       |      yang1 | (linalg = :deterministic,) |         0.48 |       5.30 |     35.33 |   3.43 |       0.28 |  0.01 |   44.83 |
|       |      yang1 |    (linalg = :randomized,) |         0.53 |       5.12 |      2.27 |   3.30 |       0.37 |  0.01 |   11.60 |
|       |      yang1 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.15 |    0.61 |
|       |   bayes148 | (linalg = :deterministic,) |         0.35 |       5.12 |     22.23 |   3.35 |       0.30 |  0.01 |   31.38 |
|       |   bayes148 |    (linalg = :randomized,) |         0.60 |       4.74 |      6.36 |   3.32 |       0.30 |  0.29 |   15.61 |
|       |   bayes148 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.64 |
|       |     mayr42 | (linalg = :deterministic,) |         0.31 |       1.72 |      7.80 |  25.85 |       0.23 |  0.01 |   35.92 |
|       |     mayr42 |    (linalg = :randomized,) |         0.31 |       2.29 |      0.71 |  24.92 |       0.23 |  0.01 |   28.46 |
|       |     mayr42 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.08 |
|       |   rand-4-4 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.00 |
|       |   rand-4-4 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.00 |
|       |   rand-4-4 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.00 |
|       |   rand-4-5 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.02 |   0.00 |       0.00 |  0.00 |    0.03 |
|       |   rand-4-5 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |    0.02 |
|       |   rand-4-5 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.01 |
|       |   rand-5-4 | (linalg = :deterministic,) |         0.00 |       0.01 |      0.15 |   0.00 |       0.00 |  0.01 |    0.18 |
|       |   rand-5-4 |    (linalg = :randomized,) |         0.00 |       0.01 |      0.05 |   0.00 |       0.00 |  0.01 |    0.08 |
|       |   rand-5-4 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.04 |
|       |   rand-5-5 | (linalg = :deterministic,) |         0.05 |       0.18 |      2.35 |   0.01 |       0.01 |  0.08 |    2.68 |
|       |   rand-5-5 |    (linalg = :randomized,) |         0.02 |       0.22 |      0.69 |   0.01 |       0.01 |  0.10 |    1.06 |
|       |   rand-5-5 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.01 |    0.63 |
:-------+------------+----------------------------+--------------+------------+-----------+--------+------------+-------+---------:
| Total |            |                            |        24.30 |      99.46 |   1648.17 | 108.40 |      28.01 | 33.80 | 2080.27 |
'-------'------------'----------------------------'--------------'------------'-----------'--------'------------'-------'---------'
                                  Count linalg ops
.-------.------------.----------------------------.-----.---------------.-----------.
|       |       Name |                    Options | % p |    addmul % p |   mul % p |
:-------+------------+----------------------------+-----+---------------+-----------:
|       |  chandra-9 | (linalg = :deterministic,) |   0 |      12929526 |     22049 |
|       |  chandra-9 |    (linalg = :randomized,) |   0 |       3793045 |     31902 |
|       |  chandra-9 |         (:groebner_apply,) |   0 |        195529 |     22049 |
|       | chandra-10 | (linalg = :deterministic,) |   0 |      81659111 |     74777 |
|       | chandra-10 |    (linalg = :randomized,) |   0 |      21626165 |    111660 |
|       | chandra-10 |         (:groebner_apply,) |   0 |        943692 |     74777 |
|       | chandra-11 | (linalg = :deterministic,) |   0 |     514499062 |    256512 |
|       | chandra-11 |    (linalg = :randomized,) |   0 |     125574185 |    411037 |
|       | chandra-11 |         (:groebner_apply,) |   0 |       4736007 |    256512 |
|       | chandra-12 | (linalg = :deterministic,) |   0 |    3237160656 |    886078 |
|       | chandra-12 |    (linalg = :randomized,) |   0 |     724197188 |   1551182 |
|       | chandra-12 |         (:groebner_apply,) |   0 |      24418740 |    886078 |
|       | chandra-13 | (linalg = :deterministic,) |   0 |   20367862985 |   3075088 |
|       | chandra-13 |    (linalg = :randomized,) |   0 |    4347498587 |   5748975 |
|       | chandra-13 |         (:groebner_apply,) |   0 |     128369221 |   3075088 |
|       |   cyclic-7 | (linalg = :deterministic,) |   0 |      72996366 |    115935 |
|       |   cyclic-7 |    (linalg = :randomized,) |   0 |      18008616 |    123495 |
|       |   cyclic-7 |         (:groebner_apply,) |   0 |      10351699 |    115935 |
|       |   cyclic-8 | (linalg = :deterministic,) |   0 |    2092236718 |   1007207 |
|       |   cyclic-8 |    (linalg = :randomized,) |   0 |     433566324 |   1081294 |
|       |   cyclic-8 |         (:groebner_apply,) |   0 |     283983858 |   1007207 |
|       |     eco-11 | (linalg = :deterministic,) |   0 |     475504869 |     70944 |
|       |     eco-11 |    (linalg = :randomized,) |   0 |      81926626 |    145169 |
|       |     eco-11 |         (:groebner_apply,) |   0 |      39072150 |     70944 |
|       |     eco-12 | (linalg = :deterministic,) |   0 |    3842762167 |    230898 |
|       |     eco-12 |    (linalg = :randomized,) |   0 |     551625580 |    518708 |
|       |     eco-12 |         (:groebner_apply,) |   0 |     309145352 |    230898 |
|       |     eco-13 | (linalg = :deterministic,) |   0 |   32573400714 |    819993 |
|       |     eco-13 |    (linalg = :randomized,) |   0 |    2867593017 |   1820466 |
|       |     eco-13 |         (:groebner_apply,) |   0 |    1358702260 |    819993 |
|       |     eco-14 | (linalg = :deterministic,) |   0 |  292334414838 |   3120702 |
|       |     eco-14 |    (linalg = :randomized,) |   0 |   31875417544 |   7151240 |
|       |     eco-14 |         (:groebner_apply,) |   0 |   20930018214 |   3120702 |
|       |     noon-7 | (linalg = :deterministic,) |   0 |       8136307 |     65987 |
|       |     noon-7 |    (linalg = :randomized,) |   0 |      13328386 |    220984 |
|       |     noon-7 |         (:groebner_apply,) |   0 |        673467 |     65987 |
|       |     noon-8 | (linalg = :deterministic,) |   0 |      71866242 |    392602 |
|       |     noon-8 |    (linalg = :randomized,) |   0 |     157654102 |   1647909 |
|       |     noon-8 |         (:groebner_apply,) |   0 |       5116538 |    392602 |
|       |     noon-9 | (linalg = :deterministic,) |   0 |     621970346 |   2344062 |
|       |     noon-9 |    (linalg = :randomized,) |   0 |    1891707215 |  12564524 |
|       |     noon-9 |         (:groebner_apply,) |   0 |      38139478 |   2344062 |
|       |    noon-10 | (linalg = :deterministic,) |   0 |    5385499271 |  14131850 |
|       |    noon-10 |    (linalg = :randomized,) |   0 |   23021232575 |  95247409 |
|       |    noon-10 |         (:groebner_apply,) |   0 |     281629370 |  14131850 |
|       |  henrion-6 | (linalg = :deterministic,) |   0 |      11917980 |     33227 |
|       |  henrion-6 |    (linalg = :randomized,) |   0 |       5117342 |     35423 |
|       |  henrion-6 |         (:groebner_apply,) |   0 |       1305119 |     33227 |
|       |  henrion-7 | (linalg = :deterministic,) |   0 |    2347580853 |   1073152 |
|       |  henrion-7 |    (linalg = :randomized,) |   0 |     608389170 |   1145067 |
|       |  henrion-7 |         (:groebner_apply,) |   0 |     230728956 |   1073152 |
|       |  henrion-8 | (linalg = :deterministic,) |   0 |  840184897086 |  49706073 |
|       |  henrion-8 |    (linalg = :randomized,) |   0 |  149262888866 |  52472559 |
|       |  henrion-8 |         (:groebner_apply,) |   0 |   77244283686 |  49706073 |
|       |  katsura-9 | (linalg = :deterministic,) |   0 |     318008502 |     98102 |
|       |  katsura-9 |    (linalg = :randomized,) |   0 |      41038011 |    102447 |
|       |  katsura-9 |         (:groebner_apply,) |   0 |      16554189 |     98102 |
|       | katsura-10 | (linalg = :deterministic,) |   0 |    2668201940 |    382812 |
|       | katsura-10 |    (linalg = :randomized,) |   0 |     279521918 |    397101 |
|       | katsura-10 |         (:groebner_apply,) |   0 |     126114961 |    382812 |
|       | katsura-11 | (linalg = :deterministic,) |   0 |   22439551011 |   1480142 |
|       | katsura-11 |    (linalg = :randomized,) |   0 |    1961456047 |   1536029 |
|       | katsura-11 |         (:groebner_apply,) |   0 |     967598341 |   1480142 |
|       | katsura-12 | (linalg = :deterministic,) |   0 |  189523146948 |   5800299 |
|       | katsura-12 |    (linalg = :randomized,) |   0 |   14229436358 |   6017690 |
|       | katsura-12 |         (:groebner_apply,) |   0 |    7565970968 |   5800299 |
|       |   reimer-6 | (linalg = :deterministic,) |   0 |      18980777 |     24766 |
|       |   reimer-6 |    (linalg = :randomized,) |   0 |       5225166 |     51903 |
|       |   reimer-6 |         (:groebner_apply,) |   0 |        780712 |     24766 |
|       |   reimer-7 | (linalg = :deterministic,) |   0 |     564039977 |    195808 |
|       |   reimer-7 |    (linalg = :randomized,) |   0 |     118079783 |    386868 |
|       |   reimer-7 |         (:groebner_apply,) |   0 |      18186252 |    195808 |
|       |   reimer-8 | (linalg = :deterministic,) |   0 |   19934877599 |   1683422 |
|       |   reimer-8 |    (linalg = :randomized,) |   0 |    3407446512 |   3368766 |
|       |   reimer-8 |         (:groebner_apply,) |   0 |     547790127 |   1683422 |
|       |    hexapod | (linalg = :deterministic,) |   0 |       1336148 |      6343 |
|       |    hexapod |    (linalg = :randomized,) |   0 |        550056 |      6884 |
|       |    hexapod |         (:groebner_apply,) |   0 |        225999 |      6343 |
|       |    cholera | (linalg = :deterministic,) |   0 |   46309780158 |   7937478 |
|       |    cholera |    (linalg = :randomized,) |   0 |   14743806559 |   9064124 |
|       |    cholera |         (:groebner_apply,) |   0 |   10008476493 |   7937478 |
|       |       hiv2 | (linalg = :deterministic,) |   0 |     256388895 |    357427 |
|       |       hiv2 |    (linalg = :randomized,) |   0 |     253835403 |    708496 |
|       |       hiv2 |         (:groebner_apply,) |   0 |      57216026 |    357427 |
|       |    goodwin | (linalg = :deterministic,) |   0 |   33076197249 |  14828836 |
|       |    goodwin |    (linalg = :randomized,) |   0 |   48027040127 |  33649314 |
|       |    goodwin |         (:groebner_apply,) |   0 |    9051946088 |  14828836 |
|       |        crn | (linalg = :deterministic,) |   0 |     274725427 |     64122 |
|       |        crn |    (linalg = :randomized,) |   0 |      39712390 |     89640 |
|       |        crn |         (:groebner_apply,) |   0 |      19352094 |     64122 |
|       |      alea6 | (linalg = :deterministic,) |   0 |     108933805 |     86902 |
|       |      alea6 |    (linalg = :randomized,) |   0 |      28145967 |     88543 |
|       |      alea6 |         (:groebner_apply,) |   0 |      17786410 |     86902 |
|       |   jason210 | (linalg = :deterministic,) |   0 |      46582152 |    284900 |
|       |   jason210 |    (linalg = :randomized,) |   0 |     197061780 |    860743 |
|       |   jason210 |         (:groebner_apply,) |   0 |       5071227 |    284900 |
|       |   gametwo2 | (linalg = :deterministic,) |   0 |   27996088632 |   4153336 |
|       |   gametwo2 |    (linalg = :randomized,) |   0 |    5852613481 |   4605771 |
|       |   gametwo2 |         (:groebner_apply,) |   0 |    3974027141 |   4153336 |
|       |      yang1 | (linalg = :deterministic,) |   0 |       6017076 |      3300 |
|       |      yang1 |    (linalg = :randomized,) |   0 |     430765304 |   2125961 |
|       |      yang1 |         (:groebner_apply,) |   0 |        108900 |      3300 |
|       |   bayes148 | (linalg = :deterministic,) |   0 |      13014323 |     26436 |
|       |   bayes148 |    (linalg = :randomized,) |   0 |     130113000 |    816599 |
|       |   bayes148 |         (:groebner_apply,) |   0 |        129533 |     26436 |
|       |     mayr42 | (linalg = :deterministic,) |   0 |       2440056 |      2154 |
|       |     mayr42 |    (linalg = :randomized,) |   0 |      23144339 |    116439 |
|       |     mayr42 |         (:groebner_apply,) |   0 |         37692 |      2154 |
|       |   rand-4-4 | (linalg = :deterministic,) |   0 |       3209538 |     12688 |
|       |   rand-4-4 |    (linalg = :randomized,) |   0 |       1535907 |     12688 |
|       |   rand-4-4 |         (:groebner_apply,) |   0 |        793502 |     12688 |
|       |   rand-4-5 | (linalg = :deterministic,) |   0 |      24810933 |     55959 |
|       |   rand-4-5 |    (linalg = :randomized,) |   0 |      10936484 |     55959 |
|       |   rand-4-5 |         (:groebner_apply,) |   0 |       6468844 |     55959 |
|       |   rand-5-4 | (linalg = :deterministic,) |   0 |     206514122 |    174504 |
|       |   rand-5-4 |    (linalg = :randomized,) |   0 |      62228052 |    174504 |
|       |   rand-5-4 |         (:groebner_apply,) |   0 |      39184241 |    174504 |
|       |   rand-5-5 | (linalg = :deterministic,) |   0 |    3223310940 |   1228002 |
|       |   rand-5-5 |    (linalg = :randomized,) |   0 |     904999332 |   1228002 |
|       |   rand-5-5 |         (:groebner_apply,) |   0 |     645597871 |   1228002 |
:-------+------------+----------------------------+-----+---------------+-----------:
| Total |            |                            |   0 | 1991974518761 | 480123222 |
'-------'------------'----------------------------'-----'---------------'-----------'

                                                        Timings in seconds
.-------.------------.----------------------------.--------------.------------.-----------.--------.------------.-------.--------.
|       |       Name |                    Options | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |    IO |  Total |
:-------+------------+----------------------------+--------------+------------+-----------+--------+------------+-------+--------:
|       |  chandra-9 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.01 |   0.01 |       0.00 |  0.00 |   0.02 |
|       |  chandra-9 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |  chandra-9 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       | chandra-10 | (linalg = :deterministic,) |         0.01 |       0.00 |      0.07 |   0.01 |       0.00 |  0.01 |   0.10 |
|       | chandra-10 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.02 |   0.01 |       0.00 |  0.01 |   0.05 |
|       | chandra-10 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.02 |   0.02 |
|       | chandra-11 | (linalg = :deterministic,) |         0.01 |       0.10 |      0.44 |   0.03 |       0.01 |  0.02 |   0.61 |
|       | chandra-11 |    (linalg = :randomized,) |         0.03 |       0.10 |      0.12 |   0.03 |       0.01 |  0.02 |   0.31 |
|       | chandra-11 |         (:groebner_apply,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.02 |   0.03 |
|       | chandra-12 | (linalg = :deterministic,) |         0.05 |       0.19 |      2.67 |   0.10 |       0.05 |  0.20 |   3.26 |
|       | chandra-12 |    (linalg = :randomized,) |         0.07 |       0.09 |      0.64 |   0.21 |       0.02 |  0.10 |   1.13 |
|       | chandra-12 |         (:groebner_apply,) |         0.00 |       0.01 |      0.04 |   0.00 |       0.02 |  0.08 |   0.15 |
|       | chandra-13 | (linalg = :deterministic,) |         0.23 |       0.58 |     17.65 |   0.43 |       0.08 |  0.66 |  19.63 |
|       | chandra-13 |    (linalg = :randomized,) |         0.27 |       0.59 |      4.28 |   0.60 |       0.10 |  0.67 |   6.50 |
|       | chandra-13 |         (:groebner_apply,) |         0.00 |       0.05 |      0.17 |   0.00 |       0.06 |  0.33 |   0.61 |
|       |   cyclic-7 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.06 |   0.01 |       0.00 |  0.00 |   0.07 |
|       |   cyclic-7 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.02 |   0.01 |       0.00 |  0.00 |   0.03 |
|       |   cyclic-7 |         (:groebner_apply,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |   cyclic-8 | (linalg = :deterministic,) |         0.03 |       0.03 |      1.62 |   0.03 |       0.00 |  0.01 |   1.72 |
|       |   cyclic-8 |    (linalg = :randomized,) |         0.03 |       0.03 |      0.37 |   0.04 |       0.00 |  0.00 |   0.48 |
|       |   cyclic-8 |         (:groebner_apply,) |         0.00 |       0.11 |      0.23 |   0.00 |       0.00 |  0.00 |   0.34 |
|       |     eco-11 | (linalg = :deterministic,) |         0.01 |       0.02 |      0.43 |   0.01 |       0.02 |  0.01 |   0.51 |
|       |     eco-11 |    (linalg = :randomized,) |         0.01 |       0.02 |      0.08 |   0.01 |       0.03 |  0.00 |   0.14 |
|       |     eco-11 |         (:groebner_apply,) |         0.00 |       0.01 |      0.12 |   0.00 |       0.02 |  0.00 |   0.15 |
|       |     eco-12 | (linalg = :deterministic,) |         0.03 |       0.18 |      3.29 |   0.02 |       0.14 |  0.03 |   3.69 |
|       |     eco-12 |    (linalg = :randomized,) |         0.12 |       0.11 |      0.44 |   0.02 |       0.14 |  0.03 |   0.85 |
|       |     eco-12 |         (:groebner_apply,) |         0.00 |       0.10 |      0.19 |   0.00 |       0.14 |  0.11 |   0.54 |
|       |     eco-13 | (linalg = :deterministic,) |         0.22 |       0.42 |     29.30 |   0.08 |       0.13 |  0.31 |  30.46 |
|       |     eco-13 |    (linalg = :randomized,) |         0.09 |       0.54 |      2.86 |   0.08 |       0.04 |  0.21 |   3.83 |
|       |     eco-13 |         (:groebner_apply,) |         0.00 |       0.34 |      1.15 |   0.00 |       0.04 |  0.18 |   1.71 |
|       |     noon-7 | (linalg = :deterministic,) |         0.00 |       0.02 |      0.03 |   0.01 |       0.00 |  0.00 |   0.06 |
|       |     noon-7 |    (linalg = :randomized,) |         0.04 |       0.02 |      0.03 |   0.01 |       0.00 |  0.00 |   0.09 |
|       |     noon-7 |         (:groebner_apply,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.02 |
|       |     noon-8 | (linalg = :deterministic,) |         0.03 |       0.10 |      0.39 |   0.06 |       0.02 |  0.03 |   0.63 |
|       |     noon-8 |    (linalg = :randomized,) |         0.03 |       0.11 |      0.26 |   0.05 |       0.02 |  0.03 |   0.50 |
|       |     noon-8 |         (:groebner_apply,) |         0.00 |       0.02 |      0.04 |   0.00 |       0.02 |  0.00 |   0.08 |
|       |     noon-9 | (linalg = :deterministic,) |         0.15 |       0.74 |      2.82 |   0.39 |       0.14 |  0.41 |   4.65 |
|       |     noon-9 |    (linalg = :randomized,) |         0.14 |       0.91 |      2.77 |   0.38 |       0.14 |  0.27 |   4.61 |
|       |     noon-9 |         (:groebner_apply,) |         0.00 |       0.13 |      0.27 |   0.00 |       0.10 |  0.24 |   0.73 |
|       |    noon-10 | (linalg = :deterministic,) |         0.97 |       6.19 |     27.32 |   3.27 |       1.54 |  1.84 |  41.12 |
|       |    noon-10 |    (linalg = :randomized,) |         0.93 |       5.74 |     31.75 |   3.37 |       1.43 |  2.02 |  45.23 |
|       |    noon-10 |         (:groebner_apply,) |         0.00 |       0.89 |      2.07 |   0.00 |       0.65 |  1.33 |   4.94 |
|       |  henrion-6 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.02 |
|       |  henrion-6 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |  henrion-6 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |  henrion-7 | (linalg = :deterministic,) |         0.03 |       0.50 |      1.83 |   0.01 |       0.01 |  0.09 |   2.48 |
|       |  henrion-7 |    (linalg = :randomized,) |         0.03 |       0.29 |      0.52 |   0.01 |       0.01 |  0.10 |   0.97 |
|       |  henrion-7 |         (:groebner_apply,) |         0.00 |       0.17 |      0.20 |   0.00 |       0.01 |  0.01 |   0.39 |
|       |  katsura-9 | (linalg = :deterministic,) |         0.01 |       0.01 |      0.27 |   0.00 |       0.00 |  0.11 |   0.39 |
|       |  katsura-9 |    (linalg = :randomized,) |         0.01 |       0.01 |      0.04 |   0.00 |       0.00 |  0.00 |   0.06 |
|       |  katsura-9 |         (:groebner_apply,) |         0.00 |       0.00 |      0.02 |   0.00 |       0.00 |  0.00 |   0.02 |
|       | katsura-10 | (linalg = :deterministic,) |         0.02 |       0.03 |      2.16 |   0.01 |       0.01 |  0.15 |   2.37 |
|       | katsura-10 |    (linalg = :randomized,) |         0.02 |       0.03 |      0.25 |   0.01 |       0.01 |  0.11 |   0.43 |
|       | katsura-10 |         (:groebner_apply,) |         0.00 |       0.02 |      0.11 |   0.00 |       0.00 |  0.02 |   0.15 |
|       | katsura-11 | (linalg = :deterministic,) |         0.07 |       0.13 |     17.89 |   0.04 |       0.12 |  0.26 |  18.50 |
|       | katsura-11 |    (linalg = :randomized,) |         0.07 |       0.24 |      1.80 |   0.04 |       0.02 |  0.28 |   2.44 |
|       | katsura-11 |         (:groebner_apply,) |         0.00 |       0.07 |      0.77 |   0.00 |       0.02 |  0.04 |   0.90 |
|       | katsura-12 | (linalg = :deterministic,) |         0.51 |       0.91 |    159.90 |   0.15 |       0.07 |  1.24 | 162.79 |
|       | katsura-12 |    (linalg = :randomized,) |         0.30 |       0.76 |     12.62 |   0.16 |       0.16 |  1.00 |  14.99 |
|       | katsura-12 |         (:groebner_apply,) |         0.00 |       0.68 |      6.02 |   0.00 |       0.07 |  0.46 |   7.22 |
|       |   reimer-6 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.02 |   0.00 |       0.00 |  0.00 |   0.03 |
|       |   reimer-6 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.02 |
|       |   reimer-6 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |   reimer-7 | (linalg = :deterministic,) |         0.01 |       0.09 |      0.48 |   0.00 |       0.00 |  0.00 |   0.59 |
|       |   reimer-7 |    (linalg = :randomized,) |         0.01 |       0.10 |      0.12 |   0.00 |       0.00 |  0.01 |   0.24 |
|       |   reimer-7 |         (:groebner_apply,) |         0.00 |       0.01 |      0.02 |   0.00 |       0.00 |  0.00 |   0.04 |
|       |   reimer-8 | (linalg = :deterministic,) |         0.23 |       2.22 |     16.92 |   0.02 |       0.03 |  0.13 |  19.55 |
|       |   reimer-8 |    (linalg = :randomized,) |         0.12 |       2.02 |      3.31 |   0.02 |       0.03 |  0.35 |   5.86 |
|       |   reimer-8 |         (:groebner_apply,) |         0.00 |       0.39 |      0.49 |   0.00 |       0.03 |  0.06 |   0.97 |
|       |    hexapod | (linalg = :deterministic,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |    hexapod |    (linalg = :randomized,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |    hexapod |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |    cholera | (linalg = :deterministic,) |         0.99 |       0.30 |     43.20 |   1.63 |       0.00 |  0.00 |  46.12 |
|       |    cholera |    (linalg = :randomized,) |         0.96 |       0.53 |     13.35 |   1.50 |       0.00 |  0.00 |  16.34 |
|       |    cholera |         (:groebner_apply,) |         0.00 |       0.33 |      8.34 |   0.00 |       0.00 |  0.00 |   8.67 |
|       |       hiv2 | (linalg = :deterministic,) |         0.05 |       0.12 |      0.43 |   0.42 |       0.00 |  0.00 |   1.02 |
|       |       hiv2 |    (linalg = :randomized,) |         0.06 |       0.12 |      0.32 |   0.34 |       0.00 |  0.00 |   0.84 |
|       |       hiv2 |         (:groebner_apply,) |         0.00 |       0.02 |      0.10 |   0.00 |       0.00 |  0.00 |   0.12 |
|       |    goodwin | (linalg = :deterministic,) |         4.25 |       2.01 |     45.10 |  14.30 |       0.00 |  0.00 |  65.67 |
|       |    goodwin |    (linalg = :randomized,) |         4.04 |       1.93 |     41.65 |  14.26 |       0.00 |  0.00 |  61.88 |
|       |    goodwin |         (:groebner_apply,) |         0.00 |       0.71 |      7.71 |   0.00 |       0.00 |  0.00 |   8.42 |
|       |        crn | (linalg = :deterministic,) |         0.01 |       0.01 |      0.30 |   0.01 |       0.00 |  0.00 |   0.34 |
|       |        crn |    (linalg = :randomized,) |         0.01 |       0.01 |      0.05 |   0.01 |       0.00 |  0.00 |   0.09 |
|       |        crn |         (:groebner_apply,) |         0.00 |       0.01 |      0.02 |   0.00 |       0.00 |  0.00 |   0.03 |
|       |      alea6 | (linalg = :deterministic,) |         0.00 |       0.01 |      0.09 |   0.00 |       0.00 |  0.00 |   0.11 |
|       |      alea6 |    (linalg = :randomized,) |         0.00 |       0.01 |      0.03 |   0.00 |       0.00 |  0.00 |   0.05 |
|       |      alea6 |         (:groebner_apply,) |         0.00 |       0.01 |      0.04 |   0.00 |       0.00 |  0.00 |   0.04 |
|       |   jason210 | (linalg = :deterministic,) |         0.09 |       0.70 |      0.83 |   0.04 |       0.18 |  0.04 |   1.89 |
|       |   jason210 |    (linalg = :randomized,) |         0.09 |       0.81 |      0.74 |   0.14 |       0.18 |  0.03 |   1.99 |
|       |   jason210 |         (:groebner_apply,) |         0.00 |       0.08 |      0.07 |   0.00 |       0.10 |  0.00 |   0.25 |
|       |   gametwo2 | (linalg = :deterministic,) |         0.08 |       0.19 |     20.86 |   0.10 |       0.38 |  0.59 |  22.21 |
|       |   gametwo2 |    (linalg = :randomized,) |         0.10 |       0.11 |      4.59 |   0.11 |       0.41 |  0.72 |   6.03 |
|       |   gametwo2 |         (:groebner_apply,) |         0.00 |       0.29 |      2.59 |   0.00 |       0.38 |  0.27 |   3.54 |
|       |      yang1 | (linalg = :deterministic,) |         0.34 |       4.66 |     35.14 |   3.17 |       0.27 |  0.14 |  43.72 |
|       |      yang1 |    (linalg = :randomized,) |         0.50 |       4.71 |      2.51 |   2.96 |       0.27 |  0.01 |  10.96 |
|       |      yang1 |         (:groebner_apply,) |         0.00 |       0.01 |      0.32 |   0.00 |       0.12 |  0.00 |   0.45 |
|       |   bayes148 | (linalg = :deterministic,) |         0.35 |       4.67 |     22.08 |   3.56 |       0.30 |  0.01 |  30.98 |
|       |   bayes148 |    (linalg = :randomized,) |         0.35 |       4.56 |      6.37 |   3.51 |       0.31 |  0.02 |  15.13 |
|       |   bayes148 |         (:groebner_apply,) |         0.00 |       0.04 |      0.12 |   0.00 |       0.19 |  0.00 |   0.35 |
|       |     mayr42 | (linalg = :deterministic,) |         0.31 |       1.55 |      7.88 |  22.54 |       0.23 |  0.01 |  32.53 |
|       |     mayr42 |    (linalg = :randomized,) |         0.31 |       1.55 |      0.71 |  22.54 |       0.23 |  0.01 |  25.34 |
|       |     mayr42 |         (:groebner_apply,) |         0.00 |       0.01 |      0.02 |   0.00 |       0.05 |  0.00 |   0.08 |
|       |   rand-4-4 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |   rand-4-4 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |   rand-4-4 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |   rand-4-5 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.02 |   0.00 |       0.00 |  0.00 |   0.03 |
|       |   rand-4-5 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.02 |
|       |   rand-4-5 |         (:groebner_apply,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |   rand-5-4 | (linalg = :deterministic,) |         0.00 |       0.01 |      0.15 |   0.00 |       0.00 |  0.01 |   0.18 |
|       |   rand-5-4 |    (linalg = :randomized,) |         0.00 |       0.01 |      0.05 |   0.00 |       0.00 |  0.01 |   0.08 |
|       |   rand-5-4 |         (:groebner_apply,) |         0.00 |       0.01 |      0.03 |   0.00 |       0.00 |  0.00 |   0.04 |
|       |   rand-5-5 | (linalg = :deterministic,) |         0.02 |       0.09 |      2.45 |   0.01 |       0.01 |  0.08 |   2.65 |
|       |   rand-5-5 |    (linalg = :randomized,) |         0.02 |       0.23 |      0.69 |   0.01 |       0.01 |  0.09 |   1.05 |
|       |   rand-5-5 |         (:groebner_apply,) |         0.00 |       0.06 |      0.48 |   0.00 |       0.01 |  0.06 |   0.61 |
:-------+------------+----------------------------+--------------+------------+-----------+--------+------------+-------+--------:
| Total |            |                            |        17.91 |      57.70 |    629.33 | 100.89 |       9.35 | 15.75 | 830.96 |
'-------'------------'----------------------------'--------------'------------'-----------'--------'------------'-------'--------'
                                                       Count linalg ops
.-------.------------.----------------------------.--------.----------------------------.---------------------.--------------.
|       |       Name |                    Options |    % p |                 addmul % p |             mul % p | rows reduced |
:-------+------------+----------------------------+--------+----------------------------+---------------------+--------------:
|       |  chandra-9 | (linalg = :deterministic,) | (0, 0) |         (309096, 12929526) |        (291, 22049) |    (2333, 0) |
|       |  chandra-9 |    (linalg = :randomized,) | (0, 0) |          (116967, 3793045) |        (291, 31902) |     (907, 0) |
|       |  chandra-9 |         (:groebner_apply,) | (0, 0) |             (4803, 195529) |        (291, 22049) |     (838, 0) |
|       | chandra-10 | (linalg = :deterministic,) | (0, 0) |        (1364071, 81659111) |        (556, 74777) |    (5028, 0) |
|       | chandra-10 |    (linalg = :randomized,) | (0, 0) |         (468861, 21626165) |       (556, 111660) |    (1727, 0) |
|       | chandra-10 |         (:groebner_apply,) | (0, 0) |            (15127, 943692) |        (556, 74777) |    (1624, 0) |
|       | chandra-11 | (linalg = :deterministic,) | (0, 0) |       (6042156, 514499062) |      (1078, 256512) |   (10872, 0) |
|       | chandra-11 |    (linalg = :randomized,) | (0, 0) |       (1938232, 125574185) |      (1078, 411037) |    (3337, 0) |
|       | chandra-11 |         (:groebner_apply,) | (0, 0) |           (49764, 4736007) |      (1078, 256512) |    (3180, 0) |
|       | chandra-12 | (linalg = :deterministic,) | (0, 0) |     (26809397, 3237160656) |      (2113, 886078) |   (23514, 0) |
|       | chandra-12 |    (linalg = :randomized,) | (0, 0) |       (7887100, 724197188) |     (2113, 1551182) |    (6508, 0) |
|       | chandra-12 |         (:groebner_apply,) | (0, 0) |         (167971, 24418740) |      (2113, 886078) |    (6274, 0) |
|       | chandra-13 | (linalg = :deterministic,) | (0, 0) |   (119293828, 20367862985) |     (4173, 3075088) |   (50763, 0) |
|       | chandra-13 |    (linalg = :randomized,) | (0, 0) |     (33534789, 4347498587) |     (4173, 5748975) |   (12797, 0) |
|       | chandra-13 |         (:groebner_apply,) | (0, 0) |        (575809, 128369221) |     (4173, 3075088) |   (12442, 0) |
|       |   cyclic-7 | (linalg = :deterministic,) | (0, 0) |         (569754, 72996366) |       (560, 115935) |    (3575, 0) |
|       |   cyclic-7 |    (linalg = :randomized,) | (0, 0) |         (158268, 18008616) |       (585, 123495) |    (1503, 0) |
|       |   cyclic-7 |         (:groebner_apply,) | (0, 0) |          (71797, 10351699) |       (560, 115935) |    (1379, 0) |
|       |   cyclic-8 | (linalg = :deterministic,) | (0, 0) |      (6368280, 2092236718) |     (1451, 1007207) |    (9641, 0) |
|       |   cyclic-8 |    (linalg = :randomized,) | (0, 0) |       (1435136, 433566324) |     (1509, 1081294) |    (3666, 0) |
|       |   cyclic-8 |         (:groebner_apply,) | (0, 0) |        (741639, 283983858) |     (1451, 1007207) |    (3426, 0) |
|       |     eco-11 | (linalg = :deterministic,) | (0, 0) |       (6813603, 475504869) |        (206, 70944) |    (6786, 0) |
|       |     eco-11 |    (linalg = :randomized,) | (0, 0) |        (1136949, 81926626) |       (477, 145169) |    (3928, 0) |
|       |     eco-11 |         (:groebner_apply,) | (0, 0) |         (429259, 39072150) |        (206, 70944) |    (3814, 0) |
|       |     eco-12 | (linalg = :deterministic,) | (0, 0) |     (37407616, 3842762167) |       (347, 230898) |   (14727, 0) |
|       |     eco-12 |    (linalg = :randomized,) | (0, 0) |       (4867500, 551625580) |       (876, 518708) |    (8603, 0) |
|       |     eco-12 |         (:groebner_apply,) | (0, 0) |       (2242028, 309145352) |       (347, 230898) |    (8437, 0) |
|       |     eco-13 | (linalg = :deterministic,) | (0, 0) |   (199651778, 32573400714) |       (625, 819993) |   (18678, 0) |
|       |     eco-13 |    (linalg = :randomized,) | (0, 0) |     (18557873, 2867593017) |     (1627, 1820466) |    (4987, 0) |
|       |     eco-13 |         (:groebner_apply,) | (0, 0) |      (8120972, 1358702260) |       (625, 819993) |    (4747, 0) |
|       |     noon-7 | (linalg = :deterministic,) | (0, 0) |          (570374, 8136307) |        (431, 65987) |    (3873, 0) |
|       |     noon-7 |    (linalg = :randomized,) | (0, 0) |         (577272, 13328386) |       (488, 220984) |    (1571, 0) |
|       |     noon-7 |         (:groebner_apply,) | (0, 0) |            (52440, 673467) |        (431, 65987) |    (1471, 0) |
|       |     noon-8 | (linalg = :deterministic,) | (0, 0) |        (4151066, 71866242) |      (1204, 392602) |   (11891, 0) |
|       |     noon-8 |    (linalg = :randomized,) | (0, 0) |       (4500473, 157654102) |     (1330, 1647909) |    (4179, 0) |
|       |     noon-8 |         (:groebner_apply,) | (0, 0) |          (329109, 5116538) |      (1204, 392602) |    (3998, 0) |
|       |     noon-9 | (linalg = :deterministic,) | (0, 0) |      (29551911, 621970346) |     (3384, 2344062) |   (36546, 0) |
|       |     noon-9 |    (linalg = :randomized,) | (0, 0) |     (34188366, 1891707215) |    (3673, 12564524) |   (11353, 0) |
|       |     noon-9 |         (:groebner_apply,) | (0, 0) |        (2024701, 38139478) |     (3384, 2344062) |   (11028, 0) |
|       |    noon-10 | (linalg = :deterministic,) | (0, 0) |    (210017594, 5385499271) |    (9578, 14131850) |  (112470, 0) |
|       |    noon-10 |    (linalg = :randomized,) | (0, 0) |   (261933493, 23021232575) |   (10263, 95247409) |   (31376, 0) |
|       |    noon-10 |         (:groebner_apply,) | (0, 0) |      (12235856, 281629370) |    (9578, 14131850) |   (30799, 0) |
|       |  henrion-6 | (linalg = :deterministic,) | (0, 0) |         (101121, 11917980) |         (89, 33227) |     (510, 0) |
|       |  henrion-6 |    (linalg = :randomized,) | (0, 0) |           (51812, 5117342) |         (89, 35423) |     (308, 0) |
|       |  henrion-6 |         (:groebner_apply,) | (0, 0) |           (12053, 1305119) |         (89, 33227) |     (268, 0) |
|       |  henrion-7 | (linalg = :deterministic,) | (0, 0) |      (5523833, 2347580853) |      (414, 1073152) |    (2845, 0) |
|       |  henrion-7 |    (linalg = :randomized,) | (0, 0) |       (1678330, 608389170) |      (414, 1145067) |    (1345, 0) |
|       |  henrion-7 |         (:groebner_apply,) | (0, 0) |        (575628, 230728956) |      (414, 1073152) |    (1243, 0) |
|       |  katsura-9 | (linalg = :deterministic,) | (0, 0) |       (3469245, 318008502) |        (264, 98102) |    (2525, 0) |
|       |  katsura-9 |    (linalg = :randomized,) | (0, 0) |         (500114, 41038011) |       (264, 102447) |     (875, 0) |
|       |  katsura-9 |         (:groebner_apply,) | (0, 0) |         (175599, 16554189) |        (264, 98102) |     (807, 0) |
|       | katsura-10 | (linalg = :deterministic,) | (0, 0) |     (17761984, 2668201940) |       (528, 382812) |    (5552, 0) |
|       | katsura-10 |    (linalg = :randomized,) | (0, 0) |       (2112967, 279521918) |       (528, 397101) |    (1707, 0) |
|       | katsura-10 |         (:groebner_apply,) | (0, 0) |        (815068, 126114961) |       (528, 382812) |    (1601, 0) |
|       | katsura-11 | (linalg = :deterministic,) | (0, 0) |    (92264634, 22439551011) |     (1040, 1480142) |   (11974, 0) |
|       | katsura-11 |    (linalg = :randomized,) | (0, 0) |      (9083763, 1961456047) |     (1040, 1536029) |    (3299, 0) |
|       | katsura-11 |         (:groebner_apply,) | (0, 0) |       (3855561, 967598341) |     (1040, 1480142) |    (3139, 0) |
|       | katsura-12 | (linalg = :deterministic,) | (0, 0) |  (459707479, 189523146948) |     (2080, 5800299) |   (25990, 0) |
|       | katsura-12 |    (linalg = :randomized,) | (0, 0) |    (39685710, 14229436358) |     (2080, 6017690) |    (6503, 0) |
|       | katsura-12 |         (:groebner_apply,) | (0, 0) |     (17955929, 7565970968) |     (2080, 5800299) |    (6261, 0) |
|       |   reimer-6 | (linalg = :deterministic,) | (0, 0) |         (183521, 18980777) |        (133, 24766) |    (1023, 0) |
|       |   reimer-6 |    (linalg = :randomized,) | (0, 0) |           (61940, 5225166) |        (198, 51903) |     (552, 0) |
|       |   reimer-6 |         (:groebner_apply,) | (0, 0) |             (9560, 780712) |        (133, 24766) |     (498, 0) |
|       |   reimer-7 | (linalg = :deterministic,) | (0, 0) |       (1997944, 564039977) |       (286, 195808) |    (2491, 0) |
|       |   reimer-7 |    (linalg = :randomized,) | (0, 0) |        (531725, 118079783) |       (422, 386868) |    (1180, 0) |
|       |   reimer-7 |         (:groebner_apply,) | (0, 0) |          (78192, 18186252) |       (286, 195808) |    (1084, 0) |
|       |   reimer-8 | (linalg = :deterministic,) | (0, 0) |    (25595958, 19934877599) |      (638, 1683422) |    (6687, 0) |
|       |   reimer-8 |    (linalg = :randomized,) | (0, 0) |      (5538241, 3407446512) |      (948, 3368766) |    (2686, 0) |
|       |   reimer-8 |         (:groebner_apply,) | (0, 0) |        (832301, 547790127) |      (638, 1683422) |    (2524, 0) |
|       |    hexapod | (linalg = :deterministic,) | (0, 0) |           (33627, 1336148) |         (111, 6343) |     (574, 0) |
|       |    hexapod |    (linalg = :randomized,) | (0, 0) |            (13248, 550056) |         (113, 6884) |     (313, 0) |
|       |    hexapod |         (:groebner_apply,) | (0, 0) |             (5218, 225999) |         (111, 6343) |     (282, 0) |
|       |    cholera | (linalg = :deterministic,) | (0, 0) |   (429632132, 46309780158) |     (7263, 7937478) |   (52050, 0) |
|       |    cholera |    (linalg = :randomized,) | (0, 0) |    (70178440, 14743806559) |     (7412, 9064124) |   (15249, 0) |
|       |    cholera |         (:groebner_apply,) | (0, 0) |    (35594937, 10008476493) |     (7263, 7937478) |   (14868, 0) |
|       |       hiv2 | (linalg = :deterministic,) | (0, 0) |       (7974891, 256388895) |      (3274, 357427) |   (26115, 0) |
|       |       hiv2 |    (linalg = :randomized,) | (0, 0) |       (5016874, 253835403) |      (3521, 708496) |    (7446, 0) |
|       |       hiv2 |         (:groebner_apply,) | (0, 0) |         (935474, 57216026) |      (3274, 357427) |    (7097, 0) |
|       |    goodwin | (linalg = :deterministic,) | (0, 0) |   (778939675, 33076197249) |   (12705, 14828836) |  (130621, 0) |
|       |    goodwin |    (linalg = :randomized,) | (0, 0) |    (86370588, 48027040127) |   (13034, 33649314) |   (26901, 0) |
|       |    goodwin |         (:groebner_apply,) | (0, 0) |     (15185403, 9051946088) |   (12705, 14828836) |   (26109, 0) |
|       |        crn | (linalg = :deterministic,) | (0, 0) |       (4998265, 274725427) |        (693, 64122) |    (5721, 0) |
|       |        crn |    (linalg = :randomized,) | (0, 0) |         (854896, 39712390) |        (776, 89640) |    (1700, 0) |
|       |        crn |         (:groebner_apply,) | (0, 0) |         (348346, 19352094) |        (693, 64122) |    (1597, 0) |
|       |      alea6 | (linalg = :deterministic,) | (0, 0) |       (1142542, 108933805) |        (227, 86902) |    (1673, 0) |
|       |      alea6 |    (linalg = :randomized,) | (0, 0) |         (307470, 28145967) |        (227, 88543) |    (1005, 0) |
|       |      alea6 |         (:groebner_apply,) | (0, 0) |         (183965, 17786410) |        (227, 86902) |     (941, 0) |
|       |   jason210 | (linalg = :deterministic,) | (0, 0) |        (2398111, 46582152) |       (863, 284900) |    (7555, 0) |
|       |   jason210 |    (linalg = :randomized,) | (0, 0) |       (9968139, 197061780) |       (897, 860743) |    (2950, 0) |
|       |   jason210 |         (:groebner_apply,) | (0, 0) |          (253964, 5071227) |       (863, 284900) |    (2694, 0) |
|       |   gametwo2 | (linalg = :deterministic,) | (0, 0) |    (22038890, 27996088632) |     (2587, 4153336) |   (16137, 0) |
|       |   gametwo2 |    (linalg = :randomized,) | (0, 0) |      (5078474, 5852613481) |     (2587, 4605771) |    (6932, 0) |
|       |   gametwo2 |         (:groebner_apply,) | (0, 0) |      (3513278, 3974027141) |     (2587, 4153336) |    (6727, 0) |
|       |      yang1 | (linalg = :deterministic,) | (0, 0) |          (787948, 6017076) |         (300, 3300) |  (103996, 0) |
|       |      yang1 |    (linalg = :randomized,) | (0, 0) |       (5469641, 430765304) |     (4695, 2125961) |   (14477, 0) |
|       |      yang1 |         (:groebner_apply,) | (0, 0) |            (17460, 108900) |         (300, 3300) |   (14151, 0) |
|       |   bayes148 | (linalg = :deterministic,) | (0, 0) |        (3473479, 13014323) |        (912, 26436) |   (64737, 0) |
|       |   bayes148 |    (linalg = :randomized,) | (0, 0) |      (26405535, 130113000) |      (3605, 816599) |   (11300, 0) |
|       |   bayes148 |         (:groebner_apply,) | (0, 0) |            (44406, 129533) |        (912, 26436) |   (10836, 0) |
|       |     mayr42 | (linalg = :deterministic,) | (0, 0) |         (1220028, 2440056) |        (2154, 2154) |  (356243, 0) |
|       |     mayr42 |    (linalg = :randomized,) | (0, 0) |        (7224564, 23144339) |      (8490, 116439) |   (26914, 0) |
|       |     mayr42 |         (:groebner_apply,) | (0, 0) |             (18846, 37692) |        (2154, 2154) |   (25514, 0) |
|       |   rand-4-4 | (linalg = :deterministic,) | (0, 0) |           (20486, 3209538) |         (59, 12688) |     (283, 0) |
|       |   rand-4-4 |    (linalg = :randomized,) | (0, 0) |           (11023, 1535907) |         (59, 12688) |     (200, 0) |
|       |   rand-4-4 |         (:groebner_apply,) | (0, 0) |             (5272, 793502) |         (59, 12688) |     (178, 0) |
|       |   rand-4-5 | (linalg = :deterministic,) | (0, 0) |          (68369, 24810933) |        (107, 55959) |     (514, 0) |
|       |   rand-4-5 |    (linalg = :randomized,) | (0, 0) |          (34473, 10936484) |        (107, 55959) |     (354, 0) |
|       |   rand-4-5 |         (:groebner_apply,) | (0, 0) |           (18443, 6468844) |        (107, 55959) |     (322, 0) |
|       |   rand-5-4 | (linalg = :deterministic,) | (0, 0) |        (385327, 206514122) |       (205, 174504) |    (1170, 0) |
|       |   rand-5-4 |    (linalg = :randomized,) | (0, 0) |         (144324, 62228052) |       (205, 174504) |     (663, 0) |
|       |   rand-5-4 |         (:groebner_apply,) | (0, 0) |          (74760, 39184241) |       (205, 174504) |     (616, 0) |
|       |   rand-5-5 | (linalg = :deterministic,) | (0, 0) |      (2105747, 3223310940) |      (472, 1228002) |    (2723, 0) |
|       |   rand-5-5 |    (linalg = :randomized,) | (0, 0) |        (748968, 904999332) |      (472, 1228002) |    (1494, 0) |
|       |   rand-5-5 |         (:groebner_apply,) | (0, 0) |        (432496, 645597871) |      (472, 1228002) |    (1417, 0) |
:-------+------------+----------------------------+--------+----------------------------+---------------------+--------------:
| Total |            |                            | (0, 0) | (3267121732, 580142598527) | (208024, 314845873) | (1597432, 0) |
'-------'------------'----------------------------'--------'----------------------------'---------------------'--------------'

                                                     Timings in seconds
.-------.------------.----------------------------.--------------.------------.-----------.--------.------------.-------.--------.
|       |       Name |                    Options | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |    IO |  Total |
:-------+------------+----------------------------+--------------+------------+-----------+--------+------------+-------+--------:
|       |  chandra-9 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.01 |   0.01 |       0.00 |  0.00 |   0.02 |
|       |  chandra-9 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |  chandra-9 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       | chandra-10 | (linalg = :deterministic,) |         0.00 |       0.09 |      0.09 |   0.01 |       0.00 |  0.00 |   0.19 |
|       | chandra-10 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.02 |   0.01 |       0.00 |  0.01 |   0.05 |
|       | chandra-10 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.09 |   0.09 |
|       | chandra-11 | (linalg = :deterministic,) |         0.01 |       0.04 |      0.43 |   0.03 |       0.01 |  0.04 |   0.55 |
|       | chandra-11 |    (linalg = :randomized,) |         0.01 |       0.02 |      0.21 |   0.03 |       0.01 |  0.04 |   0.32 |
|       | chandra-11 |         (:groebner_apply,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.01 |  0.01 |   0.03 |
|       | chandra-12 | (linalg = :deterministic,) |         0.05 |       0.10 |      2.65 |   0.10 |       0.05 |  0.30 |   3.24 |
|       | chandra-12 |    (linalg = :randomized,) |         0.05 |       0.19 |      0.66 |   0.10 |       0.04 |  0.18 |   1.21 |
|       | chandra-12 |         (:groebner_apply,) |         0.00 |       0.01 |      0.07 |   0.00 |       0.01 |  0.03 |   0.13 |
|       | chandra-13 | (linalg = :deterministic,) |         0.23 |       0.65 |     17.41 |   0.44 |       0.07 |  0.65 |  19.46 |
|       | chandra-13 |    (linalg = :randomized,) |         0.28 |       0.58 |      4.44 |   0.61 |       0.10 |  0.67 |   6.67 |
|       | chandra-13 |         (:groebner_apply,) |         0.00 |       0.05 |      0.18 |   0.00 |       0.06 |  0.34 |   0.63 |
|       |   cyclic-7 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.06 |   0.01 |       0.00 |  0.00 |   0.07 |
|       |   cyclic-7 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.02 |   0.01 |       0.00 |  0.00 |   0.03 |
|       |   cyclic-7 |         (:groebner_apply,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |   cyclic-8 | (linalg = :deterministic,) |         0.03 |       0.04 |      1.62 |   0.03 |       0.00 |  0.01 |   1.73 |
|       |   cyclic-8 |    (linalg = :randomized,) |         0.03 |       0.03 |      0.39 |   0.03 |       0.00 |  0.00 |   0.48 |
|       |   cyclic-8 |         (:groebner_apply,) |         0.00 |       0.12 |      0.23 |   0.00 |       0.00 |  0.00 |   0.34 |
|       |     eco-11 | (linalg = :deterministic,) |         0.01 |       0.02 |      0.43 |   0.01 |       0.02 |  0.01 |   0.49 |
|       |     eco-11 |    (linalg = :randomized,) |         0.01 |       0.02 |      0.08 |   0.01 |       0.02 |  0.01 |   0.14 |
|       |     eco-11 |         (:groebner_apply,) |         0.00 |       0.01 |      0.03 |   0.00 |       0.02 |  0.00 |   0.06 |
|       |     eco-12 | (linalg = :deterministic,) |         0.03 |       0.19 |      3.26 |   0.02 |       0.13 |  0.03 |   3.67 |
|       |     eco-12 |    (linalg = :randomized,) |         0.03 |       0.11 |      0.54 |   0.02 |       0.13 |  0.04 |   0.88 |
|       |     eco-12 |         (:groebner_apply,) |         0.00 |       0.11 |      0.19 |   0.00 |       0.14 |  0.00 |   0.44 |
|       |     eco-13 | (linalg = :deterministic,) |         0.10 |       0.61 |     28.85 |   0.08 |       0.14 |  0.32 |  30.10 |
|       |     eco-13 |    (linalg = :randomized,) |         0.10 |       0.52 |      2.93 |   0.08 |       0.14 |  0.22 |   3.98 |
|       |     eco-13 |         (:groebner_apply,) |         0.00 |       0.49 |      1.13 |   0.00 |       0.04 |  0.03 |   1.69 |
|       |     noon-7 | (linalg = :deterministic,) |         0.00 |       0.02 |      0.03 |   0.01 |       0.00 |  0.00 |   0.06 |
|       |     noon-7 |    (linalg = :randomized,) |         0.00 |       0.02 |      0.03 |   0.01 |       0.00 |  0.11 |   0.17 |
|       |     noon-7 |         (:groebner_apply,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |     noon-8 | (linalg = :deterministic,) |         0.03 |       0.13 |      0.28 |   0.05 |       0.02 |  0.12 |   0.62 |
|       |     noon-8 |    (linalg = :randomized,) |         0.03 |       0.09 |      0.36 |   0.05 |       0.02 |  0.02 |   0.58 |
|       |     noon-8 |         (:groebner_apply,) |         0.00 |       0.02 |      0.04 |   0.00 |       0.02 |  0.00 |   0.08 |
|       |     noon-9 | (linalg = :deterministic,) |         0.16 |       0.68 |      2.73 |   0.36 |       0.14 |  0.41 |   4.48 |
|       |     noon-9 |    (linalg = :randomized,) |         0.14 |       0.69 |      2.95 |   0.39 |       0.14 |  0.38 |   4.68 |
|       |     noon-9 |         (:groebner_apply,) |         0.00 |       0.13 |      0.27 |   0.00 |       0.10 |  0.25 |   0.74 |
|       |    noon-10 | (linalg = :deterministic,) |         1.28 |       5.77 |     27.90 |   3.27 |       1.16 |  2.18 |  41.57 |
|       |    noon-10 |    (linalg = :randomized,) |         1.28 |       5.48 |     31.72 |   3.40 |       1.37 |  2.06 |  45.31 |
|       |    noon-10 |         (:groebner_apply,) |         0.00 |       0.91 |      2.08 |   0.00 |       0.66 |  1.32 |   4.97 |
|       |  henrion-6 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.02 |
|       |  henrion-6 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |  henrion-6 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |  henrion-7 | (linalg = :deterministic,) |         0.26 |       0.29 |      1.83 |   0.01 |       0.01 |  0.08 |   2.48 |
|       |  henrion-7 |    (linalg = :randomized,) |         0.04 |       0.24 |      0.63 |   0.01 |       0.01 |  0.09 |   1.01 |
|       |  henrion-7 |         (:groebner_apply,) |         0.00 |       0.11 |      0.20 |   0.00 |       0.07 |  0.01 |   0.39 |
|       |  katsura-9 | (linalg = :deterministic,) |         0.11 |       0.01 |      0.27 |   0.00 |       0.00 |  0.00 |   0.39 |
|       |  katsura-9 |    (linalg = :randomized,) |         0.00 |       0.01 |      0.04 |   0.00 |       0.00 |  0.00 |   0.06 |
|       |  katsura-9 |         (:groebner_apply,) |         0.00 |       0.00 |      0.02 |   0.00 |       0.00 |  0.01 |   0.03 |
|       | katsura-10 | (linalg = :deterministic,) |         0.03 |       0.03 |      2.16 |   0.01 |       0.00 |  0.12 |   2.36 |
|       | katsura-10 |    (linalg = :randomized,) |         0.02 |       0.03 |      0.27 |   0.01 |       0.00 |  0.03 |   0.36 |
|       | katsura-10 |         (:groebner_apply,) |         0.00 |       0.02 |      0.11 |   0.00 |       0.00 |  0.00 |   0.13 |
|       | katsura-11 | (linalg = :deterministic,) |         0.10 |       0.23 |     18.10 |   0.04 |       0.02 |  0.30 |  18.80 |
|       | katsura-11 |    (linalg = :randomized,) |         0.07 |       0.15 |      1.75 |   0.04 |       0.02 |  0.41 |   2.43 |
|       | katsura-11 |         (:groebner_apply,) |         0.00 |       0.23 |      0.77 |   0.00 |       0.02 |  0.19 |   1.21 |
|       | katsura-12 | (linalg = :deterministic,) |         0.46 |       1.04 |    159.28 |   0.22 |       0.19 |  1.22 | 162.41 |
|       | katsura-12 |    (linalg = :randomized,) |         0.32 |       0.88 |     12.52 |   0.16 |       0.17 |  1.00 |  15.05 |
|       | katsura-12 |         (:groebner_apply,) |         0.00 |       0.67 |      5.96 |   0.00 |       0.07 |  0.46 |   7.17 |
|       |   reimer-6 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.02 |   0.00 |       0.00 |  0.00 |   0.03 |
|       |   reimer-6 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.02 |
|       |   reimer-6 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |   reimer-7 | (linalg = :deterministic,) |         0.01 |       0.11 |      0.46 |   0.00 |       0.00 |  0.01 |   0.60 |
|       |   reimer-7 |    (linalg = :randomized,) |         0.01 |       0.23 |      0.12 |   0.00 |       0.00 |  0.12 |   0.49 |
|       |   reimer-7 |         (:groebner_apply,) |         0.00 |       0.01 |      0.02 |   0.00 |       0.00 |  0.01 |   0.05 |
|       |   reimer-8 | (linalg = :deterministic,) |         0.23 |       2.18 |     16.43 |   0.02 |       0.03 |  0.22 |  19.12 |
|       |   reimer-8 |    (linalg = :randomized,) |         0.15 |       2.10 |      3.17 |   0.02 |       0.03 |  0.33 |   5.80 |
|       |   reimer-8 |         (:groebner_apply,) |         0.00 |       0.37 |      0.49 |   0.00 |       0.03 |  0.03 |   0.92 |
|       |    hexapod | (linalg = :deterministic,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |    hexapod |    (linalg = :randomized,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |    hexapod |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |    cholera | (linalg = :deterministic,) |         1.12 |       0.42 |     43.17 |   1.41 |       0.00 |  0.00 |  46.13 |
|       |    cholera |    (linalg = :randomized,) |         1.07 |       0.35 |     13.68 |   1.38 |       0.00 |  0.00 |  16.48 |
|       |    cholera |         (:groebner_apply,) |         0.00 |       0.23 |      8.48 |   0.00 |       0.00 |  0.00 |   8.71 |
|       |       hiv2 | (linalg = :deterministic,) |         0.07 |       0.12 |      0.53 |   0.30 |       0.00 |  0.00 |   1.01 |
|       |       hiv2 |    (linalg = :randomized,) |         0.06 |       0.12 |      0.33 |   0.33 |       0.00 |  0.00 |   0.84 |
|       |       hiv2 |         (:groebner_apply,) |         0.00 |       0.04 |      0.08 |   0.00 |       0.00 |  0.00 |   0.12 |
|       |    goodwin | (linalg = :deterministic,) |         4.05 |       1.93 |     45.22 |  14.25 |       0.00 |  0.00 |  65.45 |
|       |    goodwin |    (linalg = :randomized,) |         3.90 |       1.91 |     42.18 |  13.88 |       0.00 |  0.00 |  61.87 |
|       |    goodwin |         (:groebner_apply,) |         0.00 |       0.88 |      7.83 |   0.00 |       0.00 |  0.00 |   8.71 |
|       |        crn | (linalg = :deterministic,) |         0.01 |       0.01 |      0.30 |   0.01 |       0.00 |  0.00 |   0.33 |
|       |        crn |    (linalg = :randomized,) |         0.01 |       0.01 |      0.05 |   0.01 |       0.00 |  0.00 |   0.09 |
|       |        crn |         (:groebner_apply,) |         0.00 |       0.01 |      0.02 |   0.00 |       0.00 |  0.00 |   0.03 |
|       |      alea6 | (linalg = :deterministic,) |         0.00 |       0.01 |      0.10 |   0.00 |       0.00 |  0.00 |   0.12 |
|       |      alea6 |    (linalg = :randomized,) |         0.00 |       0.01 |      0.03 |   0.00 |       0.00 |  0.00 |   0.04 |
|       |      alea6 |         (:groebner_apply,) |         0.00 |       0.01 |      0.02 |   0.00 |       0.00 |  0.00 |   0.02 |
|       |   jason210 | (linalg = :deterministic,) |         0.10 |       0.73 |      0.56 |   0.04 |       0.17 |  0.01 |   1.62 |
|       |   jason210 |    (linalg = :randomized,) |         0.08 |       0.85 |      0.70 |   0.04 |       0.17 |  0.01 |   1.85 |
|       |   jason210 |         (:groebner_apply,) |         0.00 |       0.05 |      0.11 |   0.00 |       0.10 |  0.00 |   0.25 |
|       |   gametwo2 | (linalg = :deterministic,) |         0.12 |       0.21 |     20.60 |   0.10 |       0.40 |  0.45 |  21.89 |
|       |   gametwo2 |    (linalg = :randomized,) |         0.08 |       0.10 |      4.18 |   0.10 |       0.38 |  0.66 |   5.52 |
|       |   gametwo2 |         (:groebner_apply,) |         0.00 |       0.33 |      2.59 |   0.00 |       0.42 |  0.23 |   3.57 |
|       |      yang1 | (linalg = :deterministic,) |         0.34 |       4.79 |     35.60 |   3.36 |       0.27 |  0.01 |  44.37 |
|       |      yang1 |    (linalg = :randomized,) |         0.26 |       4.88 |      2.16 |   3.13 |       0.27 |  0.01 |  10.72 |
|       |      yang1 |         (:groebner_apply,) |         0.00 |       0.01 |      0.33 |   0.00 |       0.12 |  0.00 |   0.47 |
|       |   bayes148 | (linalg = :deterministic,) |         0.38 |       4.99 |     22.38 |   3.06 |       0.30 |  0.01 |  31.12 |
|       |   bayes148 |    (linalg = :randomized,) |         0.80 |       4.47 |      6.38 |   3.13 |       0.31 |  0.02 |  15.10 |
|       |   bayes148 |         (:groebner_apply,) |         0.00 |       0.04 |      0.12 |   0.00 |       0.19 |  0.00 |   0.35 |
|       |     mayr42 | (linalg = :deterministic,) |         0.32 |       1.54 |      8.28 |  23.24 |       0.23 |  0.01 |  33.62 |
|       |     mayr42 |    (linalg = :randomized,) |         0.32 |       2.07 |      0.70 |  22.50 |       0.23 |  0.01 |  25.83 |
|       |     mayr42 |         (:groebner_apply,) |         0.00 |       0.01 |      0.02 |   0.00 |       0.05 |  0.00 |   0.08 |
|       |   rand-4-4 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |   rand-4-4 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |   rand-4-4 |         (:groebner_apply,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
|       |   rand-4-5 | (linalg = :deterministic,) |         0.00 |       0.00 |      0.02 |   0.00 |       0.00 |  0.00 |   0.03 |
|       |   rand-4-5 |    (linalg = :randomized,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.02 |
|       |   rand-4-5 |         (:groebner_apply,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.01 |
|       |   rand-5-4 | (linalg = :deterministic,) |         0.00 |       0.01 |      0.15 |   0.00 |       0.00 |  0.01 |   0.18 |
|       |   rand-5-4 |    (linalg = :randomized,) |         0.00 |       0.01 |      0.05 |   0.00 |       0.00 |  0.01 |   0.08 |
|       |   rand-5-4 |         (:groebner_apply,) |         0.00 |       0.01 |      0.03 |   0.00 |       0.00 |  0.00 |   0.04 |
|       |   rand-5-5 | (linalg = :deterministic,) |         0.02 |       0.44 |      2.34 |   0.01 |       0.01 |  0.04 |   2.86 |
|       |   rand-5-5 |    (linalg = :randomized,) |         0.02 |       0.18 |      0.70 |   0.01 |       0.01 |  0.11 |   1.03 |
|       |   rand-5-5 |         (:groebner_apply,) |         0.00 |       0.06 |      0.48 |   0.00 |       0.06 |  0.01 |   0.61 |
:-------+------------+----------------------------+--------------+------------+-----------+--------+------------+-------+--------:
| Total |            |                            |        18.88 |      58.77 |    629.50 | 100.00 |       9.21 | 16.15 | 832.54 |
'-------'------------'----------------------------'--------------'------------'-----------'--------'------------'-------'--------'
                                                               Count linalg ops
.-------.------------.----------------------------.--------.----------------------------.---------------------.--------------.---------------.
|       |       Name |                    Options |    % p |                     addmul |                 mul | rows reduced | rows reducers |
:-------+------------+----------------------------+--------+----------------------------+---------------------+--------------+---------------:
|       |  chandra-9 | (linalg = :deterministic,) | (0, 0) |         (309096, 12929526) |        (291, 22049) |    (2333, 0) |          8232 |
|       |  chandra-9 |    (linalg = :randomized,) | (0, 0) |          (116967, 3793045) |        (291, 31902) |     (907, 0) |          8232 |
|       |  chandra-9 |         (:groebner_apply,) | (0, 0) |             (4803, 195529) |        (291, 22049) |     (838, 0) |           992 |
|       | chandra-10 | (linalg = :deterministic,) | (0, 0) |        (1364071, 81659111) |        (556, 74777) |    (5028, 0) |         22173 |
|       | chandra-10 |    (linalg = :randomized,) | (0, 0) |         (468861, 21626165) |       (556, 111660) |    (1727, 0) |         22173 |
|       | chandra-10 |         (:groebner_apply,) | (0, 0) |            (15127, 943692) |        (556, 74777) |    (1624, 0) |          2169 |
|       | chandra-11 | (linalg = :deterministic,) | (0, 0) |       (6042156, 514499062) |      (1078, 256512) |   (10872, 0) |         61550 |
|       | chandra-11 |    (linalg = :randomized,) | (0, 0) |       (1938232, 125574185) |      (1078, 411037) |    (3337, 0) |         61550 |
|       | chandra-11 |         (:groebner_apply,) | (0, 0) |           (49764, 4736007) |      (1078, 256512) |    (3180, 0) |          4757 |
|       | chandra-12 | (linalg = :deterministic,) | (0, 0) |     (26809397, 3237160656) |      (2113, 886078) |   (23514, 0) |        173722 |
|       | chandra-12 |    (linalg = :randomized,) | (0, 0) |       (7887100, 724197188) |     (2113, 1551182) |    (6508, 0) |        173722 |
|       | chandra-12 |         (:groebner_apply,) | (0, 0) |         (167971, 24418740) |      (2113, 886078) |    (6274, 0) |         10420 |
|       | chandra-13 | (linalg = :deterministic,) | (0, 0) |   (119293828, 20367862985) |     (4173, 3075088) |   (50763, 0) |        499786 |
|       | chandra-13 |    (linalg = :randomized,) | (0, 0) |     (33534789, 4347498587) |     (4173, 5748975) |   (12797, 0) |        499786 |
|       | chandra-13 |         (:groebner_apply,) | (0, 0) |        (575809, 128369221) |     (4173, 3075088) |   (12442, 0) |         22742 |
|       |   cyclic-7 | (linalg = :deterministic,) | (0, 0) |         (569754, 72996366) |       (560, 115935) |    (3575, 0) |          8926 |
|       |   cyclic-7 |    (linalg = :randomized,) | (0, 0) |         (158268, 18008616) |       (585, 123495) |    (1503, 0) |          8926 |
|       |   cyclic-7 |         (:groebner_apply,) | (0, 0) |          (71797, 10351699) |       (560, 115935) |    (1379, 0) |          4018 |
|       |   cyclic-8 | (linalg = :deterministic,) | (0, 0) |      (6368280, 2092236718) |     (1451, 1007207) |    (9641, 0) |         53051 |
|       |   cyclic-8 |    (linalg = :randomized,) | (0, 0) |       (1435136, 433566324) |     (1509, 1081294) |    (3666, 0) |         53051 |
|       |   cyclic-8 |         (:groebner_apply,) | (0, 0) |        (741639, 283983858) |     (1451, 1007207) |    (3426, 0) |         20616 |
|       |     eco-11 | (linalg = :deterministic,) | (0, 0) |       (6813603, 475504869) |        (206, 70944) |    (6786, 0) |         68548 |
|       |     eco-11 |    (linalg = :randomized,) | (0, 0) |        (1136949, 81926626) |       (477, 145169) |    (3928, 0) |         68548 |
|       |     eco-11 |         (:groebner_apply,) | (0, 0) |         (429259, 39072150) |        (206, 70944) |    (3814, 0) |         33342 |
|       |     eco-12 | (linalg = :deterministic,) | (0, 0) |     (37407616, 3842762167) |       (347, 230898) |   (14727, 0) |        188457 |
|       |     eco-12 |    (linalg = :randomized,) | (0, 0) |       (4867500, 551625580) |       (876, 518708) |    (8603, 0) |        188457 |
|       |     eco-12 |         (:groebner_apply,) | (0, 0) |       (2242028, 309145352) |       (347, 230898) |    (8437, 0) |        124132 |
|       |     eco-13 | (linalg = :deterministic,) | (0, 0) |   (199651778, 32573400714) |       (625, 819993) |   (18678, 0) |        544624 |
|       |     eco-13 |    (linalg = :randomized,) | (0, 0) |     (18557873, 2867593017) |     (1627, 1820466) |    (4987, 0) |        544624 |
|       |     eco-13 |         (:groebner_apply,) | (0, 0) |      (8120972, 1358702260) |       (625, 819993) |    (4747, 0) |        369746 |
|       |     noon-7 | (linalg = :deterministic,) | (0, 0) |          (570374, 8136307) |        (431, 65987) |    (3873, 0) |         78113 |
|       |     noon-7 |    (linalg = :randomized,) | (0, 0) |         (577272, 13328386) |       (488, 220984) |    (1571, 0) |         78113 |
|       |     noon-7 |         (:groebner_apply,) | (0, 0) |            (52440, 673467) |        (431, 65987) |    (1471, 0) |         19831 |
|       |     noon-8 | (linalg = :deterministic,) | (0, 0) |        (4151066, 71866242) |      (1204, 392602) |   (11891, 0) |        382756 |
|       |     noon-8 |    (linalg = :randomized,) | (0, 0) |       (4500473, 157654102) |     (1330, 1647909) |    (4179, 0) |        382756 |
|       |     noon-8 |         (:groebner_apply,) | (0, 0) |          (329109, 5116538) |      (1204, 392602) |    (3998, 0) |        100149 |
|       |     noon-9 | (linalg = :deterministic,) | (0, 0) |      (29551911, 621970346) |     (3384, 2344062) |   (36546, 0) |       1833901 |
|       |     noon-9 |    (linalg = :randomized,) | (0, 0) |     (34188366, 1891707215) |    (3673, 12564524) |   (11353, 0) |       1833901 |
|       |     noon-9 |         (:groebner_apply,) | (0, 0) |        (2024701, 38139478) |     (3384, 2344062) |   (11028, 0) |        491872 |
|       |    noon-10 | (linalg = :deterministic,) | (0, 0) |    (210017594, 5385499271) |    (9578, 14131850) |  (112470, 0) |       8664489 |
|       |    noon-10 |    (linalg = :randomized,) | (0, 0) |   (261933493, 23021232575) |   (10263, 95247409) |   (31376, 0) |       8664489 |
|       |    noon-10 |         (:groebner_apply,) | (0, 0) |      (12235856, 281629370) |    (9578, 14131850) |   (30799, 0) |       2372557 |
|       |  henrion-6 | (linalg = :deterministic,) | (0, 0) |         (101121, 11917980) |         (89, 33227) |     (510, 0) |          9156 |
|       |  henrion-6 |    (linalg = :randomized,) | (0, 0) |           (51812, 5117342) |         (89, 35423) |     (308, 0) |          9156 |
|       |  henrion-6 |         (:groebner_apply,) | (0, 0) |           (12053, 1305119) |         (89, 33227) |     (268, 0) |          4200 |
|       |  henrion-7 | (linalg = :deterministic,) | (0, 0) |      (5523833, 2347580853) |      (414, 1073152) |    (2845, 0) |        120673 |
|       |  henrion-7 |    (linalg = :randomized,) | (0, 0) |       (1678330, 608389170) |      (414, 1145067) |    (1345, 0) |        120673 |
|       |  henrion-7 |         (:groebner_apply,) | (0, 0) |        (575628, 230728956) |      (414, 1073152) |    (1243, 0) |         72504 |
|       |  katsura-9 | (linalg = :deterministic,) | (0, 0) |       (3469245, 318008502) |        (264, 98102) |    (2525, 0) |         25572 |
|       |  katsura-9 |    (linalg = :randomized,) | (0, 0) |         (500114, 41038011) |       (264, 102447) |     (875, 0) |         25572 |
|       |  katsura-9 |         (:groebner_apply,) | (0, 0) |         (175599, 16554189) |        (264, 98102) |     (807, 0) |         12931 |
|       | katsura-10 | (linalg = :deterministic,) | (0, 0) |     (17761984, 2668201940) |       (528, 382812) |    (5552, 0) |         69492 |
|       | katsura-10 |    (linalg = :randomized,) | (0, 0) |       (2112967, 279521918) |       (528, 397101) |    (1707, 0) |         69492 |
|       | katsura-10 |         (:groebner_apply,) | (0, 0) |        (815068, 126114961) |       (528, 382812) |    (1601, 0) |         32822 |
|       | katsura-11 | (linalg = :deterministic,) | (0, 0) |    (92264634, 22439551011) |     (1040, 1480142) |   (11974, 0) |        179164 |
|       | katsura-11 |    (linalg = :randomized,) | (0, 0) |      (9083763, 1961456047) |     (1040, 1536029) |    (3299, 0) |        179164 |
|       | katsura-11 |         (:groebner_apply,) | (0, 0) |       (3855561, 967598341) |     (1040, 1480142) |    (3139, 0) |         96847 |
|       | katsura-12 | (linalg = :deterministic,) | (0, 0) |  (459707479, 189523146948) |     (2080, 5800299) |   (25990, 0) |        469984 |
|       | katsura-12 |    (linalg = :randomized,) | (0, 0) |    (39685710, 14229436358) |     (2080, 6017690) |    (6503, 0) |        469984 |
|       | katsura-12 |         (:groebner_apply,) | (0, 0) |     (17955929, 7565970968) |     (2080, 5800299) |    (6261, 0) |        255037 |
|       |   reimer-6 | (linalg = :deterministic,) | (0, 0) |         (183521, 18980777) |        (133, 24766) |    (1023, 0) |          9495 |
|       |   reimer-6 |    (linalg = :randomized,) | (0, 0) |           (61940, 5225166) |        (198, 51903) |     (552, 0) |          9495 |
|       |   reimer-6 |         (:groebner_apply,) | (0, 0) |             (9560, 780712) |        (133, 24766) |     (498, 0) |          2369 |
|       |   reimer-7 | (linalg = :deterministic,) | (0, 0) |       (1997944, 564039977) |       (286, 195808) |    (2491, 0) |         62998 |
|       |   reimer-7 |    (linalg = :randomized,) | (0, 0) |        (531725, 118079783) |       (422, 386868) |    (1180, 0) |         62998 |
|       |   reimer-7 |         (:groebner_apply,) | (0, 0) |          (78192, 18186252) |       (286, 195808) |    (1084, 0) |         13353 |
|       |   reimer-8 | (linalg = :deterministic,) | (0, 0) |    (25595958, 19934877599) |      (638, 1683422) |    (6687, 0) |        451940 |
|       |   reimer-8 |    (linalg = :randomized,) | (0, 0) |      (5538241, 3407446512) |      (948, 3368766) |    (2686, 0) |        451940 |
|       |   reimer-8 |         (:groebner_apply,) | (0, 0) |        (832301, 547790127) |      (638, 1683422) |    (2524, 0) |         74652 |
|       |    hexapod | (linalg = :deterministic,) | (0, 0) |           (33627, 1336148) |         (111, 6343) |     (574, 0) |          1581 |
|       |    hexapod |    (linalg = :randomized,) | (0, 0) |            (13248, 550056) |         (113, 6884) |     (313, 0) |          1581 |
|       |    hexapod |         (:groebner_apply,) | (0, 0) |             (5218, 225999) |         (111, 6343) |     (282, 0) |           461 |
|       |    cholera | (linalg = :deterministic,) | (0, 0) |   (429632132, 46309780158) |     (7263, 7937478) |   (52050, 0) |        297513 |
|       |    cholera |    (linalg = :randomized,) | (0, 0) |    (70178440, 14743806559) |     (7412, 9064124) |   (15249, 0) |        297513 |
|       |    cholera |         (:groebner_apply,) | (0, 0) |    (35594937, 10008476493) |     (7263, 7937478) |   (14868, 0) |         78988 |
|       |       hiv2 | (linalg = :deterministic,) | (0, 0) |       (7974891, 256388895) |      (3274, 357427) |   (26115, 0) |        133457 |
|       |       hiv2 |    (linalg = :randomized,) | (0, 0) |       (5016874, 253835403) |      (3521, 708496) |    (7446, 0) |        133457 |
|       |       hiv2 |         (:groebner_apply,) | (0, 0) |         (935474, 57216026) |      (3274, 357427) |    (7097, 0) |         29061 |
|       |    goodwin | (linalg = :deterministic,) | (0, 0) |   (778939675, 33076197249) |   (12705, 14828836) |  (130621, 0) |       1912786 |
|       |    goodwin |    (linalg = :randomized,) | (0, 0) |    (86370588, 48027040127) |   (13034, 33649314) |   (26901, 0) |       1912786 |
|       |    goodwin |         (:groebner_apply,) | (0, 0) |     (15185403, 9051946088) |   (12705, 14828836) |   (26109, 0) |        115108 |
|       |        crn | (linalg = :deterministic,) | (0, 0) |       (4998265, 274725427) |        (693, 64122) |    (5721, 0) |         17029 |
|       |        crn |    (linalg = :randomized,) | (0, 0) |         (854896, 39712390) |        (776, 89640) |    (1700, 0) |         17029 |
|       |        crn |         (:groebner_apply,) | (0, 0) |         (348346, 19352094) |        (693, 64122) |    (1597, 0) |          7233 |
|       |      alea6 | (linalg = :deterministic,) | (0, 0) |       (1142542, 108933805) |        (227, 86902) |    (1673, 0) |         22978 |
|       |      alea6 |    (linalg = :randomized,) | (0, 0) |         (307470, 28145967) |        (227, 88543) |    (1005, 0) |         22978 |
|       |      alea6 |         (:groebner_apply,) | (0, 0) |         (183965, 17786410) |        (227, 86902) |     (941, 0) |         13339 |
|       |   jason210 | (linalg = :deterministic,) | (0, 0) |        (2398111, 46582152) |       (863, 284900) |    (7555, 0) |       1237652 |
|       |   jason210 |    (linalg = :randomized,) | (0, 0) |       (9968139, 197061780) |       (897, 860743) |    (2950, 0) |       1237652 |
|       |   jason210 |         (:groebner_apply,) | (0, 0) |          (253964, 5071227) |       (863, 284900) |    (2694, 0) |        184150 |
|       |   gametwo2 | (linalg = :deterministic,) | (0, 0) |    (22038890, 27996088632) |     (2587, 4153336) |   (16137, 0) |         34643 |
|       |   gametwo2 |    (linalg = :randomized,) | (0, 0) |      (5078474, 5852613481) |     (2587, 4605771) |    (6932, 0) |         34643 |
|       |   gametwo2 |         (:groebner_apply,) | (0, 0) |      (3513278, 3974027141) |     (2587, 4153336) |    (6727, 0) |         22635 |
|       |      yang1 | (linalg = :deterministic,) | (0, 0) |          (787948, 6017076) |         (300, 3300) |  (103996, 0) |        540323 |
|       |      yang1 |    (linalg = :randomized,) | (0, 0) |       (5469641, 430765304) |     (4695, 2125961) |   (14477, 0) |        540323 |
|       |      yang1 |         (:groebner_apply,) | (0, 0) |            (17460, 108900) |         (300, 3300) |   (14151, 0) |         17460 |
|       |   bayes148 | (linalg = :deterministic,) | (0, 0) |        (3473479, 13014323) |        (912, 26436) |   (64737, 0) |       3564965 |
|       |   bayes148 |    (linalg = :randomized,) | (0, 0) |      (26405535, 130113000) |      (3605, 816599) |   (11300, 0) |       3564965 |
|       |   bayes148 |         (:groebner_apply,) | (0, 0) |            (44406, 129533) |        (912, 26436) |   (10836, 0) |         44079 |
|       |     mayr42 | (linalg = :deterministic,) | (0, 0) |         (1220028, 2440056) |        (2154, 2154) |  (356243, 0) |        773203 |
|       |     mayr42 |    (linalg = :randomized,) | (0, 0) |        (7224564, 23144339) |      (8490, 116439) |   (26914, 0) |        773203 |
|       |     mayr42 |         (:groebner_apply,) | (0, 0) |             (18846, 37692) |        (2154, 2154) |   (25514, 0) |         16156 |
|       |   rand-4-4 | (linalg = :deterministic,) | (0, 0) |           (20486, 3209538) |         (59, 12688) |     (283, 0) |          1850 |
|       |   rand-4-4 |    (linalg = :randomized,) | (0, 0) |           (11023, 1535907) |         (59, 12688) |     (200, 0) |          1850 |
|       |   rand-4-4 |         (:groebner_apply,) | (0, 0) |             (5272, 793502) |         (59, 12688) |     (178, 0) |           986 |
|       |   rand-4-5 | (linalg = :deterministic,) | (0, 0) |          (68369, 24810933) |        (107, 55959) |     (514, 0) |          3967 |
|       |   rand-4-5 |    (linalg = :randomized,) | (0, 0) |          (34473, 10936484) |        (107, 55959) |     (354, 0) |          3967 |
|       |   rand-4-5 |         (:groebner_apply,) | (0, 0) |           (18443, 6468844) |        (107, 55959) |     (322, 0) |          2439 |
|       |   rand-5-4 | (linalg = :deterministic,) | (0, 0) |        (385327, 206514122) |       (205, 174504) |    (1170, 0) |          9352 |
|       |   rand-5-4 |    (linalg = :randomized,) | (0, 0) |         (144324, 62228052) |       (205, 174504) |     (663, 0) |          9352 |
|       |   rand-5-4 |         (:groebner_apply,) | (0, 0) |          (74760, 39184241) |       (205, 174504) |     (616, 0) |          5611 |
|       |   rand-5-5 | (linalg = :deterministic,) | (0, 0) |      (2105747, 3223310940) |      (472, 1228002) |    (2723, 0) |         26404 |
|       |   rand-5-5 |    (linalg = :randomized,) | (0, 0) |        (748968, 904999332) |      (472, 1228002) |    (1494, 0) |         26404 |
|       |   rand-5-5 |         (:groebner_apply,) | (0, 0) |        (432496, 645597871) |      (472, 1228002) |    (1417, 0) |         17672 |
:-------+------------+----------------------------+--------+----------------------------+---------------------+--------------+---------------:
| Total |            |                            | (0, 0) | (3267121732, 580142598527) | (208024, 314845873) | (1597432, 0) |      49826446 |
'-------'------------'----------------------------'--------'----------------------------'---------------------'--------------'---------------'

=#
