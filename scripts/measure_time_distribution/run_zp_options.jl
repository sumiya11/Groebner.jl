using AbstractAlgebra, Groebner, BenchmarkTools, PrettyTables, Printf, TimerOutputs

include("gizmos.jl")

k = AbstractAlgebra.GF(2^30+3)
system_solving = [
    ("chandra-13", Groebner.Examples.chandran(13, k=k)),
    ("reimer-8", Groebner.Examples.reimern(8, k=k)),
    ("eco-13", Groebner.Examples.econ(13, k=k)),
    ("noon-10", Groebner.Examples.noonn(10, k=k)),
    ("katsura-12", Groebner.Examples.katsuran(12, k=k)),
]
sian = [
    ("cholera", Groebner.Examples.Cholera(k=k)),
    ("goodwin", Groebner.Examples.Goodwin_with_weights(k=k)),
]
other = [
    ("jason210", Groebner.Examples.jason210(k=k)),
    ("yang1", Groebner.Examples.yang1(k=k)),
    ("bayes148", Groebner.Examples.bayes148(k=k)),
    ("mayr42", Groebner.Examples.mayr42(k=k)),
]
random = []
systems = vcat(system_solving, sian, other, random)

println("Running the following systems: ", map(first, systems))

timers = []
for (name, sys) in systems
    @info "Running $name.."

    for kws in [
        # (_use_divmask=false, monoms=:dense),
        # (_use_divmask=true, monoms=:dense),
        # (_use_divmask=false, monoms=:packed),
        (_use_divmask=true, monoms=:packed),
    ]
        @info "Options:" kws
        for t in 1:2
            Groebner._DIVMASK_CHECKS[] = 0
            Groebner._DIVMASK_FALSE_POSITIVES[] = 0
	    Groebner._MAX_DEG[] = 0
            TimerOutputs.enable_timer!(Groebner._TIMER); reset_timer!(Groebner._TIMER);
            groebner(sys; kws...);
            if t == 2 show(Groebner._TIMER, allocations=false); println() end
        end
        push!(timers, [name, kws, copy(Groebner._TIMER), Groebner._DIVMASK_CHECKS[], Groebner._DIVMASK_FALSE_POSITIVES[], Groebner._MAX_DEG[], Groebner._DIVMASK_FALSE_POSITIVES[] / max(Groebner._DIVMASK_CHECKS[], 1)])
    end    
end

labels = ["Name", "Options", "Select Pairs", "Symb Prepr", "Reduction", "Update", "Autoreduce", "IO", "Total"]
labels2 = ["Name", "Options", "Checks", "False Positives", "Max total degree", "False Positive Rate"]
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

threshold = 0.15
hl_select = TextHighlighter((v,i,j) -> (total = v[i,end]; j == findfirst(==("Select Pairs"), labels) && v[i,j] isa Number && v[i,j] / total > threshold), crayon"red bold")
hl_symb_prepr = TextHighlighter((v,i,j) -> (total = v[i,end]; j == findfirst(==("Symb Prepr"), labels) && v[i,j] isa Number && v[i,j] / total > threshold), crayon"red bold")
hl_update = TextHighlighter((v,i,j) -> (total = v[i,end]; j == findfirst(==("Update"), labels) && v[i,j] isa Number && v[i,j] / total > threshold), crayon"red bold")
hl_io = TextHighlighter((v,i,j) -> (total = v[i,end]; j == findfirst(==("IO"), labels) && v[i,j] isa Number && v[i,j] / total > threshold), crayon"bold bg:light_gray")

pretty_table(
    matrix, 
    column_labels=labels,
    title="Timings in seconds",
    formatters=[
        (v,i,j) -> v isa Number ? @sprintf("%.2f", v / 1e9) : v, 
    ],
    # row_group_labels = [1 => "System solving", length(system_solving)+1 => "SIAN", length(system_solving)+length(sian)+1 => "Other"],
    highlighters  = [hl_select, hl_symb_prepr, hl_update, hl_io],
    table_format = TextTableFormat(borders = text_table_borders__ascii_rounded),
    summary_row_labels = ["Total"],
    summary_rows = [(data, i) -> i > 2 ? @sprintf("%.2f", sum(data[:, i]) / 1e9) : ""],
    fit_table_in_display_vertically = false,
)

matrix2 = permutedims(reduce(hcat, data2))

pretty_table(
    matrix2, 
    column_labels=labels2,
    title="Divmask Stats",
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
.-------.------------.------------------------------------------.--------------.------------.-----------.--------.------------.-------.---------.
|       |       Name |                                  Options | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |    IO |   Total |
:-------+------------+------------------------------------------+--------------+------------+-----------+--------+------------+-------+---------:
|       | chandra-13 |  (_use_divmask = false, monoms = :dense) |         1.94 |      11.78 |      8.07 |   2.07 |       0.99 |  1.64 |   26.49 |
|       | chandra-13 |   (_use_divmask = true, monoms = :dense) |         1.98 |       3.81 |      8.18 |   1.77 |       0.43 |  1.63 |   17.81 |
|       | chandra-13 | (_use_divmask = false, monoms = :packed) |         0.50 |       2.72 |      7.67 |   0.90 |       0.26 |  1.23 |   13.28 |
|       | chandra-13 |  (_use_divmask = true, monoms = :packed) |         0.65 |       1.17 |      7.76 |   0.88 |       0.18 |  1.06 |   11.71 |
|       |   reimer-8 |  (_use_divmask = false, monoms = :dense) |         1.33 |      22.01 |      5.84 |   0.06 |       0.29 |  0.58 |   30.11 |
|       |   reimer-8 |   (_use_divmask = true, monoms = :dense) |         1.16 |      20.07 |      5.93 |   0.06 |       0.18 |  1.24 |   28.66 |
|       |   reimer-8 | (_use_divmask = false, monoms = :packed) |         0.25 |       5.03 |      5.67 |   0.04 |       0.08 |  0.69 |   11.75 |
|       |   reimer-8 |  (_use_divmask = true, monoms = :packed) |         0.26 |       4.84 |      5.65 |   0.04 |       0.06 |  0.60 |   11.45 |
|       |     eco-13 |  (_use_divmask = false, monoms = :dense) |         0.88 |       4.69 |      5.35 |   0.37 |       0.21 |  0.74 |   12.24 |
|       |     eco-13 |   (_use_divmask = true, monoms = :dense) |         0.85 |       3.10 |      5.45 |   0.24 |       0.12 |  0.73 |   10.49 |
|       |     eco-13 | (_use_divmask = false, monoms = :packed) |         0.19 |       1.22 |      4.86 |   0.33 |       0.09 |  0.40 |    7.08 |
|       |     eco-13 |  (_use_divmask = true, monoms = :packed) |         0.24 |       0.95 |      4.91 |   0.13 |       0.08 |  0.45 |    6.76 |
|       |    noon-10 |  (_use_divmask = false, monoms = :dense) |        10.05 |     154.59 |     79.12 |  14.64 |      20.14 |  8.97 |  287.53 |
|       |    noon-10 |   (_use_divmask = true, monoms = :dense) |         9.68 |      26.46 |     81.84 |  12.17 |       4.20 |  8.40 |  142.78 |
|       |    noon-10 | (_use_divmask = false, monoms = :packed) |         2.52 |      62.63 |     65.73 |   7.39 |       9.31 |  4.43 |  152.01 |
|       |    noon-10 |  (_use_divmask = true, monoms = :packed) |         2.13 |      10.69 |     67.12 |   6.85 |       3.19 |  4.67 |   94.65 |
|       | katsura-12 |  (_use_divmask = false, monoms = :dense) |         4.75 |       8.51 |     26.04 |   0.60 |       0.64 |  2.88 |   43.42 |
|       | katsura-12 |   (_use_divmask = true, monoms = :dense) |         4.80 |       5.87 |     25.54 |   0.51 |       0.67 |  3.03 |   40.40 |
|       | katsura-12 | (_use_divmask = false, monoms = :packed) |         0.84 |       2.04 |     23.79 |   0.27 |       0.27 |  2.13 |   29.33 |
|       | katsura-12 |  (_use_divmask = true, monoms = :packed) |         0.85 |       1.65 |     23.97 |   0.27 |       0.17 |  1.94 |   28.83 |
|       |    cholera |  (_use_divmask = false, monoms = :dense) |         2.57 |       1.39 |     22.70 |   2.68 |       0.00 |  0.00 |   29.33 |
|       |    cholera |   (_use_divmask = true, monoms = :dense) |         2.61 |       0.45 |     22.98 |   2.47 |       0.00 |  0.00 |   28.51 |
|       |    cholera | (_use_divmask = false, monoms = :packed) |         2.53 |       1.43 |     22.88 |   2.56 |       0.00 |  0.00 |   29.41 |
|       |    cholera |  (_use_divmask = true, monoms = :packed) |         2.54 |       0.36 |     22.95 |   2.63 |       0.00 |  0.00 |   28.48 |
|       |   jason210 |  (_use_divmask = false, monoms = :dense) |         0.29 |      11.79 |      1.87 |   0.18 |       3.94 |  0.28 |   18.35 |
|       |   jason210 |   (_use_divmask = true, monoms = :dense) |         0.30 |       3.01 |      1.99 |   0.11 |       0.45 |  0.27 |    6.13 |
|       |   jason210 | (_use_divmask = false, monoms = :packed) |         0.20 |       2.87 |      1.19 |   0.08 |       0.91 |  0.02 |    5.28 |
|       |   jason210 |  (_use_divmask = true, monoms = :packed) |         0.15 |       1.35 |      1.41 |   0.07 |       0.30 |  0.07 |    3.35 |
|       |      yang1 |  (_use_divmask = false, monoms = :dense) |         0.76 |     242.48 |      4.13 |   6.58 |       8.36 |  0.04 |  262.37 |
|       |      yang1 |   (_use_divmask = true, monoms = :dense) |         0.84 |       3.36 |      4.53 |   5.31 |       0.38 |  0.03 |   14.45 |
|       |      yang1 | (_use_divmask = false, monoms = :packed) |         0.77 |     228.53 |      5.17 |   5.80 |       7.92 |  0.03 |  248.23 |
|       |      yang1 |  (_use_divmask = true, monoms = :packed) |         1.21 |       2.85 |      4.41 |   5.19 |       0.38 |  0.02 |   14.05 |
|       |   bayes148 |  (_use_divmask = false, monoms = :dense) |         0.60 |     236.10 |     13.56 |   7.34 |      15.20 |  0.05 |  272.87 |
|       |   bayes148 |   (_use_divmask = true, monoms = :dense) |         0.59 |       8.66 |     13.78 |   5.46 |       0.52 |  0.05 |   29.07 |
|       |   bayes148 | (_use_divmask = false, monoms = :packed) |         0.59 |     241.82 |     14.16 |   5.84 |      15.62 |  0.85 |  278.90 |
|       |   bayes148 |  (_use_divmask = true, monoms = :packed) |         0.58 |       9.94 |     13.34 |   5.21 |       0.51 |  0.05 |   29.64 |
:-------+------------+------------------------------------------+--------------+------------+-----------+--------+------------+-------+---------:
| Total |            |                                          |        62.98 |    1350.21 |    639.50 | 107.08 |      96.05 | 49.21 | 2305.21 |
'-------'------------'------------------------------------------'--------------'------------'-----------'--------'------------'-------'---------'


W = 8 (divmask size)

                                                              Timings in seconds
.-------.------------.-----------------------------------------.--------------.------------.-----------.--------.------------.-------.--------.
|       |       Name |                                 Options | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |    IO |  Total |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+-------+--------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |         0.75 |       2.42 |      7.63 |   0.93 |       0.26 |  1.16 |  13.15 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |         0.23 |       5.65 |      5.54 |   0.04 |       0.12 |  0.66 |  12.24 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |         0.16 |       1.26 |      4.92 |   0.13 |       0.09 |  0.62 |   7.18 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |         2.47 |      16.15 |     64.93 |   6.84 |       3.93 |  4.56 |  98.88 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |         0.79 |       2.27 |     23.55 |   0.26 |       0.19 |  2.43 |  29.49 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |         2.62 |       0.53 |     23.19 |   2.56 |       0.00 |  0.00 |  28.89 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |         9.05 |       4.28 |     73.63 |  23.26 |       0.00 |  0.00 | 110.23 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |         0.09 |       2.85 |      1.60 |   0.06 |       0.79 |  0.28 |   5.67 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |         0.85 |      33.02 |      5.13 |   5.02 |       0.85 |  0.02 |  44.91 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |         0.47 |      60.98 |     12.93 |   4.88 |       3.48 |  0.43 |  83.19 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |         0.96 |      16.66 |      1.33 |  33.27 |       1.20 |  0.01 |  53.43 |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+-------+--------:
| Total |            |                                         |        18.45 |     146.07 |    224.36 |  77.25 |      10.92 | 10.18 | 487.25 |
'-------'------------'-----------------------------------------'--------------'------------'-----------'--------'------------'-------'--------'
                                                              Divmask Stats
.-------.------------.-----------------------------------------.-------------.-----------------.------------------.---------------------.
|       |       Name |                                 Options |      Checks | False Positives | Max total degree | False Positive Rate |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |   439515631 |       155005821 |               14 |            0.352674 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |   115316321 |       112842282 |               18 |            0.978546 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |   114752334 |        48999342 |               10 |            0.427001 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |  5649053323 |       438165422 |               20 |           0.0775644 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |   108177975 |        57303048 |               14 |            0.529711 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |    38243858 |         4434710 |                7 |            0.115959 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |   418564675 |        35426535 |               19 |           0.0846381 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |   734713871 |       293628256 |               55 |             0.39965 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |  7594960775 |       927513164 |                8 |            0.122122 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |  6612231176 |      1305935273 |               18 |            0.197503 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |  1373451253 |       342434611 |               33 |            0.249324 |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
| Total |            |                                         | 23198981192 |      3721688464 |              216 |             3.53469 |
'-------'------------'-----------------------------------------'-------------'-----------------'------------------'---------------------'

W = 16 (divmask size)
                                                              Timings in seconds
.-------.------------.-----------------------------------------.--------------.------------.-----------.--------.------------.------.--------.
|       |       Name |                                 Options | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |   IO |  Total |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+------+--------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |         0.45 |       2.16 |      7.43 |   1.46 |       0.24 | 1.45 |  13.19 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |         0.22 |       5.05 |      5.47 |   0.04 |       0.09 | 0.97 |  11.83 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |         0.17 |       1.12 |      4.53 |   0.13 |       0.08 | 0.52 |   6.54 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |         2.50 |      15.72 |     63.81 |   6.66 |       2.73 | 4.66 |  96.08 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |         0.56 |       1.79 |     23.09 |   0.26 |       0.36 | 1.89 |  27.94 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |         2.21 |       0.45 |     22.56 |   2.51 |       0.00 | 0.00 |  27.73 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |         9.01 |       2.64 |     73.40 |  22.01 |       0.00 | 0.00 | 107.06 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |         0.10 |       1.74 |      1.58 |   0.07 |       0.47 | 0.28 |   4.24 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |         1.01 |       6.67 |      4.10 |   4.65 |       0.44 | 0.03 |  16.90 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |         0.49 |      12.32 |     12.70 |   5.19 |       0.73 | 0.05 |  31.48 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |         0.73 |       4.64 |      1.25 |  31.10 |       0.62 | 0.01 |  38.36 |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+------+--------:
| Total |            |                                         |        17.46 |      54.29 |    219.91 |  74.06 |       5.76 | 9.84 | 381.35 |
'-------'------------'-----------------------------------------'--------------'------------'-----------'--------'------------'------'--------'
                                                              Divmask Stats
.-------.------------.-----------------------------------------.-------------.-----------------.------------------.---------------------.
|       |       Name |                                 Options |      Checks | False Positives | Max total degree | False Positive Rate |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |   439515631 |       102974522 |               14 |            0.234291 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |   115316321 |        51340359 |               18 |            0.445213 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |   114752334 |         3631308 |               10 |           0.0316447 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |  5649053323 |       362324272 |               20 |           0.0641389 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |   108177975 |        28847587 |               14 |            0.266668 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |    38243858 |         1501772 |                7 |           0.0392683 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |   418564675 |         8283633 |               19 |           0.0197906 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |   734713871 |        70793358 |               55 |            0.096355 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |  7594960775 |        54220533 |                8 |          0.00713901 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |  6612231176 |        93860375 |               18 |            0.014195 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |  1373451253 |        81313339 |               33 |           0.0592037 |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
| Total |            |                                         | 23198981192 |       859091058 |              216 |             1.27791 |
'-------'------------'-----------------------------------------'-------------'-----------------'------------------'---------------------'

W = 32 (divmask size)
                                                              Timings in seconds
.-------.------------.-----------------------------------------.--------------.------------.-----------.--------.------------.------.--------.
|       |       Name |                                 Options | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |   IO |  Total |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+------+--------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |         0.43 |       1.25 |      8.69 |   0.96 |       0.27 | 1.18 |  12.78 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |         0.23 |       4.81 |      5.88 |   0.04 |       0.07 | 0.65 |  11.69 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |         0.18 |       1.16 |      5.22 |   0.14 |       0.27 | 0.36 |   7.33 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |         2.81 |      13.90 |     68.84 |   7.38 |       2.42 | 4.29 |  99.65 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |         0.73 |       1.57 |     24.75 |   0.26 |       0.41 | 1.84 |  29.55 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |         2.48 |       0.36 |     24.58 |   2.54 |       0.00 | 0.00 |  29.97 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |        10.18 |       2.03 |     74.62 |  21.94 |       0.00 | 0.00 | 108.77 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |         0.13 |       1.33 |      1.21 |   0.07 |       0.35 | 0.23 |   3.32 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |         0.91 |       4.16 |      4.08 |   5.19 |       0.66 | 0.02 |  15.04 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |         0.62 |       9.07 |     12.13 |   6.80 |       0.57 | 0.04 |  29.25 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |         0.50 |       1.41 |      1.26 |  30.73 |       0.40 | 0.01 |  34.32 |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+------+--------:
| Total |            |                                         |        19.20 |      41.06 |    231.25 |  76.06 |       5.42 | 8.63 | 381.65 |
'-------'------------'-----------------------------------------'--------------'------------'-----------'--------'------------'------'--------'
                                                              Divmask Stats
.-------.------------.-----------------------------------------.-------------.-----------------.------------------.---------------------.
|       |       Name |                                 Options |      Checks | False Positives | Max total degree | False Positive Rate |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |   439515631 |         1554992 |               14 |          0.00353797 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |   115316321 |         6091687 |               18 |           0.0528259 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |   114752334 |          428500 |               10 |          0.00373413 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |  5649053323 |         6742704 |               20 |           0.0011936 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |   108177975 |          251360 |               14 |          0.00232358 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |    38243858 |          419784 |                7 |           0.0109765 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |   418564675 |         1909332 |               19 |          0.00456162 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |   734713871 |        17177345 |               55 |           0.0233796 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |  7594960775 |         2303632 |                8 |         0.000303311 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |  6612231176 |          743483 |               18 |         0.000112441 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |  1373451253 |        12298334 |               33 |          0.00895433 |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
| Total |            |                                         | 23198981192 |        49921153 |              216 |            0.111903 |
'-------'------------'-----------------------------------------'-------------'-----------------'------------------'---------------------'


W = 64 (divmask size)
                                                              Timings in seconds
.-------.------------.-----------------------------------------.--------------.------------.-----------.--------.------------.------.--------.
|       |       Name |                                 Options | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |   IO |  Total |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+------+--------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |         0.80 |       1.29 |      7.72 |   0.92 |       0.39 | 1.33 |  12.46 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |         0.37 |       4.86 |      5.50 |   0.04 |       0.06 | 0.69 |  11.52 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |         0.18 |       1.09 |      4.76 |   0.13 |       0.08 | 0.52 |   6.75 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |         2.34 |      11.22 |     64.91 |   6.54 |       2.85 | 4.90 |  92.76 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |         0.88 |       1.38 |     23.54 |   0.26 |       0.16 | 1.92 |  28.14 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |         2.41 |       0.36 |     22.80 |   2.52 |       0.00 | 0.00 |  28.08 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |         8.48 |       3.06 |     75.23 |  21.12 |       0.00 | 0.00 | 107.90 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |         0.15 |       1.11 |      1.23 |   0.07 |       0.30 | 0.03 |   2.87 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |         0.65 |       4.92 |      4.32 |   5.23 |       0.44 | 0.27 |  15.83 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |         0.55 |       9.66 |     13.49 |   5.25 |       0.62 | 0.04 |  29.63 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |         1.49 |       1.23 |      2.34 |  28.16 |       0.40 | 0.01 |  33.63 |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+------+--------:
| Total |            |                                         |        18.30 |      40.16 |    225.84 |  70.21 |       5.31 | 9.71 | 369.57 |
'-------'------------'-----------------------------------------'--------------'------------'-----------'--------'------------'------'--------'
                                                              Divmask Stats
.-------.------------.-----------------------------------------.-------------.-----------------.------------------.---------------------.
|       |       Name |                                 Options |      Checks | False Positives | Max total degree | False Positive Rate |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |   439515631 |          110232 |               14 |         0.000250803 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |   115316321 |           22637 |               18 |         0.000196304 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |   114752334 |           23129 |               10 |         0.000201556 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |  5649053323 |         1761955 |               20 |         0.000311903 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |   108177975 |           40300 |               14 |         0.000372534 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |    38243858 |           93231 |                7 |           0.0024378 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |   418564675 |          487998 |               19 |          0.00116588 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |   734713871 |          275021 |               55 |         0.000374324 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |  7594960775 |               0 |                8 |                 0.0 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |  6612231176 |           28032 |               18 |          4.23942e-6 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |  1373451253 |         4164593 |               33 |          0.00303221 |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
| Total |            |                                         | 23198981192 |         7007128 |              216 |          0.00834756 |
'-------'------------'-----------------------------------------'-------------'-----------------'------------------'---------------------'

W = 128 (divmask size)
                                                              Timings in seconds
.-------.------------.-----------------------------------------.--------------.------------.-----------.--------.------------.------.--------.
|       |       Name |                                 Options | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |   IO |  Total |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+------+--------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |         0.61 |       1.65 |      7.37 |   1.03 |       0.23 | 1.06 |  11.95 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |         0.43 |       4.71 |      5.70 |   0.06 |       0.12 | 0.64 |  11.65 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |         0.21 |       1.11 |      4.51 |   0.20 |       0.08 | 0.53 |   6.63 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |         2.55 |      15.07 |     64.10 |   7.14 |       2.85 | 4.59 |  96.30 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |         0.95 |       1.50 |     23.27 |   0.30 |       0.16 | 2.07 |  28.25 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |         2.21 |       0.37 |     22.45 |   2.31 |       0.00 | 0.00 |  27.34 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |         8.88 |       2.83 |     72.56 |  21.15 |       0.00 | 0.00 | 105.43 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |         0.44 |       2.38 |      1.32 |   0.29 |       0.40 | 0.02 |   4.87 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |         0.65 |       6.60 |      3.94 |   5.74 |       0.51 | 0.37 |  17.82 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |         1.95 |      12.66 |     12.29 |   6.44 |       0.74 | 0.04 |  34.15 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |         0.60 |       1.87 |      1.21 |  33.60 |       0.43 | 0.01 |  37.73 |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+------+--------:
| Total |            |                                         |        19.48 |      50.74 |    218.74 |  78.26 |       5.53 | 9.33 | 382.12 |
'-------'------------'-----------------------------------------'--------------'------------'-----------'--------'------------'------'--------'
                                                              Divmask Stats
.-------.------------.-----------------------------------------.-------------.-----------------.------------------.---------------------.
|       |       Name |                                 Options |      Checks | False Positives | Max total degree | False Positive Rate |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |   439515631 |             123 |               14 |          2.79854e-7 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |   115316321 |               0 |               18 |                 0.0 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |   114752334 |               0 |               10 |                 0.0 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |  5649053323 |           28036 |               20 |          4.96296e-6 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |   108177975 |             134 |               14 |           1.2387e-6 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |    38243858 |            1195 |                7 |          3.12468e-5 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |   418564675 |          325723 |               19 |          0.00077819 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |   734713871 |             590 |               55 |          8.03034e-7 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |  7594960775 |               0 |                8 |                 0.0 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |  6612231176 |               0 |               18 |                 0.0 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |  1373451253 |          716459 |               33 |         0.000521649 |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
| Total |            |                                         | 23198981192 |         1072260 |              216 |          0.00133837 |
'-------'------------'-----------------------------------------'-------------'-----------------'------------------'---------------------'

strategy = :first_variables, W = 32
                                                              Timings in seconds
.-------.------------.-----------------------------------------.--------------.------------.-----------.--------.------------.------.--------.
|       |       Name |                                 Options | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |   IO |  Total |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+------+--------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |         0.81 |       1.27 |      7.77 |   1.11 |       0.17 | 1.16 |  12.29 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |         0.23 |       4.90 |      5.73 |   0.07 |       0.06 | 0.63 |  11.63 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |         0.23 |       1.12 |      4.84 |   0.13 |       0.07 | 0.54 |   6.94 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |         2.65 |      12.19 |     66.46 |   6.18 |       2.22 | 4.31 |  94.02 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |         0.77 |       1.69 |     24.00 |   0.26 |       0.17 | 2.07 |  28.96 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |         2.26 |       0.86 |     23.03 |   2.61 |       0.00 | 0.00 |  28.76 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |         9.40 |       3.75 |     74.50 |  23.70 |       0.00 | 0.00 | 111.36 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |         0.13 |       1.56 |      1.24 |   0.07 |       0.35 | 0.02 |   3.38 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |         0.53 |       8.46 |      4.13 |   5.33 |       0.47 | 0.03 |  18.96 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |         1.23 |       7.54 |     13.06 |   6.25 |       0.55 | 0.40 |  29.05 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |         0.54 |       3.01 |      2.67 |  33.98 |       0.50 | 0.01 |  40.72 |
:-------+------------+-----------------------------------------+--------------+------------+-----------+--------+------------+------+--------:
| Total |            |                                         |        18.79 |      46.34 |    227.43 |  79.68 |       4.56 | 9.19 | 386.04 |
'-------'------------'-----------------------------------------'--------------'------------'-----------'--------'------------'------'--------'
                                                              Divmask Stats
.-------.------------.-----------------------------------------.-------------.-----------------.------------------.---------------------.
|       |       Name |                                 Options |      Checks | False Positives | Max total degree | False Positive Rate |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
|       | chandra-13 | (_use_divmask = true, monoms = :packed) |   439515631 |         1554992 |               14 |          0.00353797 |
|       |   reimer-8 | (_use_divmask = true, monoms = :packed) |   115316321 |         6091687 |               18 |           0.0528259 |
|       |     eco-13 | (_use_divmask = true, monoms = :packed) |   114752334 |          428500 |               10 |          0.00373413 |
|       |    noon-10 | (_use_divmask = true, monoms = :packed) |  5649053323 |         6742704 |               20 |           0.0011936 |
|       | katsura-12 | (_use_divmask = true, monoms = :packed) |   108177975 |          251360 |               14 |          0.00232358 |
|       |    cholera | (_use_divmask = true, monoms = :packed) |    38243858 |         3947795 |                7 |            0.103227 |
|       |    goodwin | (_use_divmask = true, monoms = :packed) |   418564675 |        20128880 |               19 |           0.0480903 |
|       |   jason210 | (_use_divmask = true, monoms = :packed) |   734713871 |        17177345 |               55 |           0.0233796 |
|       |      yang1 | (_use_divmask = true, monoms = :packed) |  7594960775 |        57924952 |                8 |          0.00762676 |
|       |   bayes148 | (_use_divmask = true, monoms = :packed) |  6612231176 |          743483 |               18 |         0.000112441 |
|       |     mayr42 | (_use_divmask = true, monoms = :packed) |  1373451253 |        37621057 |               33 |           0.0273916 |
:-------+------------+-----------------------------------------+-------------+-----------------+------------------+---------------------:
| Total |            |                                         | 23198981192 |       152612755 |              216 |            0.273443 |
'-------'------------'-----------------------------------------'-------------'-----------------'------------------'---------------------'


=#
