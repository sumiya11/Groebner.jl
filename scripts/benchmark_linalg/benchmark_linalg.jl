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
    ("eco-14", Groebner.Examples.econ(14, k=k)),
    ("noon-7", Groebner.Examples.noonn(7, k=k)),
    ("noon-8", Groebner.Examples.noonn(8, k=k)),
    ("noon-9", Groebner.Examples.noonn(9, k=k)),
    ("noon-10", Groebner.Examples.noonn(10, k=k)),
    ("henrion-6", Groebner.Examples.henrion6(k=k)),
    ("henrion-7", Groebner.Examples.henrion7(k=k)),
    ("henrion-8", Groebner.Examples.henrion8(k=k)),
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
            Groebner._LINALG_MOD_P[] = Groebner._LINALG_ADDMUL_MOD_P[] = Groebner._LINALG_MUL_MOD_P[] = 0
            TimerOutputs.enable_timer!(Groebner._TIMER); reset_timer!(Groebner._TIMER);
            if kws == 2 
              groebner(sys; linalg=:randomized);
              opt = (linalg=:randomized,)
      elseif kws == 1
              groebner(sys; linalg=:deterministic);
              opt = (linalg=:deterministic,)
      elseif kws == 3
              trace, _ = groebner_learn(sys);
              Groebner._LINALG_MOD_P[] = Groebner._LINALG_ADDMUL_MOD_P[] = Groebner._LINALG_MUL_MOD_P[] = 0
              TimerOutputs.enable_timer!(Groebner._TIMER); reset_timer!(Groebner._TIMER);
              groebner_apply!(trace, sys);
              opt = (:groebner_apply,)
      else
        error("beda")
      end
            if t == 2 show(Groebner._TIMER, allocations=false); println() end
        end
        push!(timers, [name, opt, copy(Groebner._TIMER), Groebner._LINALG_MOD_P[], Groebner._LINALG_ADDMUL_MOD_P[], Groebner._LINALG_MUL_MOD_P[]])
    end    
end

Base.getindex(::Nothing,key::String) = nothing
Base.getindex(timer::TimerOutput,key::String) = haskey(timer, key) ? timer.inner_timers[key] : nothing
TimerOutputs.time(::Nothing) = 0

labels = ["Name", "Options", "Select Pairs", "Symb Prepr", "Reduction", "Update", "Autoreduce", "IO", "Total"]
labels2 = ["Name", "Options", "% p", "addmul % p", "mul % p"]
data = []
data2 = []
for entry in timers
    name = entry[1]
    kws = entry[2]
    timer = entry[3]
    stats = entry[4:end]
    time_select_pairs = TimerOutputs.time(timer["f4!"]["f4_select_critical_pairs!"])
    time_symb_preprc = TimerOutputs.time(timer["f4!"]["f4_symbolic_preprocessing!"])
    time_reduction = TimerOutputs.time(timer["f4!"]["f4_reduction!"])
    time_update = TimerOutputs.time(timer["f4!"]["f4_update!"])
    time_autoreduce = TimerOutputs.time(timer["f4!"]["f4_autoreduce!"])
    time_io = TimerOutputs.time(timer["io_convert_polynomials_to_ir"]) +
              TimerOutputs.time(timer["io_convert_ir_to_polynomials"]) +
              TimerOutputs.time(timer["ir_convert_internal_to_ir"]) +
                TimerOutputs.time(timer["ir_convert_ir_to_internal"])
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
    summary_rows = [(data, i) -> i > 2 ? sum(data[:, i]) : ""],
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


=#
