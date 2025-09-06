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
    ("cyclic-9", Groebner.Examples.cyclicn(9, k=k)),
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
    for n in 4:6, d in 4:6
]), by=first)
random = sort(reduce(vcat, [
    ("rand-$(n)-$d", randsys(n, d))
    for n in 8:14, d in 2:2
]), by=first)
systems = vcat(system_solving, sian, other, random)
systems = vcat(random)

println("Running the following systems: ", map(first, systems))

timers = []
for (name, sys) in systems
    @info "Running $name.."
    # t1, t2, t3, t4, t5 = 0, 0, 0, 0, 0
    for t in 1:2
        TimerOutputs.enable_timer!(Groebner._TIMER); 
        reset_timer!(Groebner._TIMER); 
        groebner(sys); 
        if t == 2 show(Groebner._TIMER, allocations=false); println() end
    end
    push!(timers, [name, copy(Groebner._TIMER)])
end

labels = ["Name", "Select Pairs", "Symb Prepr", "Reduction", "Update", "Autoreduce", "IO", "Total"]
data = []
for entry in timers
    name = entry[1]
    timer = entry[2]
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
    push!(data, [name, time_select_pairs, time_symb_preprc, time_reduction, time_update, time_autoreduce, time_io, time_total])
end

matrix = permutedims(reduce(hcat, data))

threshold = 0.15
hl_select = TextHighlighter((v,i,j) -> (total = v[i,end]; j == 2 && v[i,j] isa Number && v[i,j] / total > threshold), crayon"red bold")
hl_symb_prepr = TextHighlighter((v,i,j) -> (total = v[i,end]; j == 3 && v[i,j] isa Number && v[i,j] / total > threshold), crayon"red bold")
hl_update = TextHighlighter((v,i,j) -> (total = v[i,end]; j == 5 && v[i,j] isa Number && v[i,j] / total > threshold), crayon"red bold")
hl_io = TextHighlighter((v,i,j) -> (total = v[i,end]; j == 7 && v[i,j] isa Number && v[i,j] / total > threshold), crayon"bold bg:light_gray")

pretty_table(
    matrix, 
    column_labels=labels,
    title="Timings in seconds",
    formatters=[
        (v,i,j) -> v isa Number ? @sprintf("%.2f", v / 1e9) : v, 
    ],
    row_group_labels = [1 => "System solving", length(system_solving)+1 => "SIAN", length(system_solving)+length(sian)+1 => "Other"],
    highlighters  = [hl_select, hl_symb_prepr, hl_update, hl_io],
    table_format = TextTableFormat(borders = text_table_borders__ascii_rounded),
    summary_row_labels = ["Total"],
    summary_rows = [(data, i) -> i > 1 ? @sprintf("%.2f", sum(data[:, i]) / 1e9) : ""],
    fit_table_in_display_vertically = false,
)

#=
                                         Timings in seconds
.-------.------------.--------------.------------.-----------.--------.------------.-------.--------.
|       |       Name | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |    IO |  Total |
:-------'------------'--------------'------------'-----------'--------'------------'-------'--------:
| System solving                                                                                    |
:-------.------------.--------------.------------.-----------.--------.------------.-------.--------:
|       |  chandra-9 |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.02 |
|       | chandra-10 |         0.01 |       0.01 |      0.04 |   0.01 |       0.00 |  0.01 |   0.09 |
|       | chandra-11 |         0.03 |       0.13 |      0.23 |   0.05 |       0.01 |  0.25 |   0.70 |
|       | chandra-12 |         0.32 |       0.21 |      1.38 |   0.22 |       0.05 |  0.51 |   2.68 |
|       | chandra-13 |         0.69 |       1.42 |      7.98 |   0.93 |       0.21 |  1.18 |  12.42 |
|       |   cyclic-7 |         0.01 |       0.00 |      0.03 |   0.01 |       0.00 |  0.00 |   0.05 |
|       |   cyclic-8 |         0.06 |       0.06 |      0.59 |   0.06 |       0.01 |  0.01 |   0.78 |
|       |   cyclic-9 |         2.13 |       4.19 |     64.07 |   1.29 |       0.09 |  0.22 |  71.99 |
|       |     eco-11 |         0.02 |       0.04 |      0.11 |   0.01 |       0.04 |  0.01 |   0.24 |
|       |     eco-12 |         0.05 |       0.36 |      0.69 |   0.04 |       0.25 |  0.07 |   1.46 |
|       |     eco-13 |         0.18 |       1.05 |      4.87 |   0.13 |       0.08 |  0.51 |   6.83 |
|       |     eco-14 |         1.04 |       4.92 |     43.99 |   0.55 |      18.34 |  1.98 |  70.83 |
|       |     noon-7 |         0.01 |       0.03 |      0.05 |   0.01 |       0.01 |  0.01 |   0.12 |
|       |     noon-8 |         0.05 |       0.20 |      0.47 |   0.09 |       0.04 |  0.06 |   0.91 |
|       |     noon-9 |         0.39 |       1.73 |      5.51 |   0.73 |       0.26 |  0.67 |   9.28 |
|       |    noon-10 |         1.98 |      11.64 |     66.48 |   7.08 |       2.82 |  4.11 |  94.11 |
|       |  henrion-6 |         0.00 |       0.01 |      0.01 |   0.00 |       0.00 |  0.00 |   0.03 |
|       |  henrion-7 |         0.27 |       0.44 |      0.87 |   0.01 |       0.04 |  0.15 |   1.78 |
|       |  henrion-8 |         5.72 |      52.64 |    225.27 |   0.35 |       1.93 | 18.37 | 304.29 |
|       |  katsura-9 |         0.01 |       0.02 |      0.06 |   0.00 |       0.00 |  0.01 |   0.11 |
|       | katsura-10 |         0.03 |       0.06 |      0.39 |   0.02 |       0.01 |  0.05 |   0.56 |
|       | katsura-11 |         0.13 |       0.26 |      2.70 |   0.06 |       0.04 |  0.16 |   3.36 |
|       | katsura-12 |         0.88 |       1.40 |     24.57 |   0.27 |       0.17 |  1.86 |  29.15 |
|       |   reimer-6 |         0.00 |       0.01 |      0.01 |   0.00 |       0.00 |  0.00 |   0.03 |
|       |   reimer-7 |         0.02 |       0.16 |      0.50 |   0.01 |       0.02 |  0.01 |   0.72 |
|       |   reimer-8 |         0.25 |       4.95 |      6.10 |   0.04 |       0.06 |  0.28 |  11.69 |
|       |    hexapod |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |   0.00 |
:-------'------------'--------------'------------'-----------'--------'------------'-------'--------:
| SIAN                                                                                              |
:-------.------------.--------------.------------.-----------.--------.------------.-------.--------:
|       |    cholera |         2.42 |       0.36 |     24.01 |   2.31 |       0.00 |  0.00 |  29.10 |
|       |       hiv2 |         0.13 |       0.08 |      0.54 |   0.42 |       0.00 |  0.00 |   1.18 |
|       |    goodwin |         9.17 |       2.45 |     75.52 |  22.47 |       0.00 |  0.00 | 109.61 |
|       |        crn |         0.02 |       0.01 |      0.07 |   0.02 |       0.00 |  0.00 |   0.13 |
:-------'------------'--------------'------------'-----------'--------'------------'-------'--------:
| Other                                                                                             |
:-------.------------.--------------.------------.-----------.--------.------------.-------.--------:
|       |      alea6 |         0.01 |       0.01 |      0.05 |   0.00 |       0.00 |  0.00 |   0.08 |
|       |   jason210 |         0.21 |       1.49 |      1.22 |   0.07 |       0.30 |  0.02 |   3.32 |
|       |   gametwo2 |         0.17 |       0.43 |      6.50 |   0.19 |       0.62 |  1.13 |   9.05 |
|       |      yang1 |         0.81 |       3.39 |      4.16 |   4.66 |       0.37 |  0.03 |  13.42 |
|       |   bayes148 |         0.58 |       7.42 |     13.07 |   6.29 |       0.52 |  0.05 |  27.94 |
|       |     mayr42 |         1.84 |       1.35 |      1.28 |  30.18 |       0.42 |  0.01 |  35.10 |
:-------'------------'--------------'------------'-----------'--------'------------'-------'--------:
| Random                                                                                             |
:-------'------------'--------------'------------'-----------'--------'------------'-------'--------:
|       |  rand-10-2 |         0.02 |       0.04 |      0.34 |   0.01 |       0.01 |  0.04 |   0.46 |
|       |  rand-11-2 |         0.09 |       0.14 |      2.05 |   0.03 |       0.03 |  0.17 |   2.50 |
|       |  rand-12-2 |         0.84 |       0.84 |     16.36 |   0.12 |       0.11 |  1.82 |  20.08 |
|       |  rand-13-2 |         3.08 |       3.29 |    103.43 |   0.43 |       0.44 |  5.65 | 116.32 |
|       |  rand-14-2 |        14.60 |      16.92 |    728.34 |   1.86 |       2.37 | 23.43 | 787.52 |
|       |   rand-8-2 |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 |  0.00 |   0.02 |
|       |   rand-9-2 |         0.01 |       0.01 |      0.06 |   0.00 |       0.00 |  0.01 |   0.09 |
|       |   rand-4-4 |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 |  0.00 |    0.01 |
|       |   rand-4-5 |         0.00 |       0.01 |      0.02 |   0.00 |       0.00 |  0.00 |    0.03 |
|       |   rand-4-6 |         0.01 |       0.03 |      0.07 |   0.00 |       0.00 |  0.02 |    0.13 |
|       |   rand-5-4 |         0.01 |       0.24 |      0.08 |   0.00 |       0.00 |  0.01 |    0.34 |
|       |   rand-5-5 |         0.04 |       0.40 |      1.04 |   0.01 |       0.02 |  0.13 |    1.64 |
|       |   rand-5-6 |         0.39 |       1.65 |     10.69 |   0.04 |       0.15 |  1.67 |   14.60 |
|       |   rand-6-4 |         0.14 |       0.45 |      3.47 |   0.02 |       0.05 |  0.52 |    4.66 |
|       |   rand-6-5 |         1.99 |       8.89 |    108.35 |   0.20 |       1.00 |  7.72 |  128.16 |
|       |   rand-6-6 |        18.12 |      92.25 |   1820.47 |   1.42 |       7.13 | 67.17 | 2006.55 |
:-------+------------+--------------+------------+-----------+--------+------------+-------+--------:
=#
