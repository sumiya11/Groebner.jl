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
    ("cholera", Groebner.Examples.Cholera(k=k)),
]
other = [
    ("jason210", Groebner.Examples.jason210(k=k)),
    ("yang1", Groebner.Examples.yang1(k=k)),
    ("bayes148", Groebner.Examples.bayes148(k=k)),
]
random = []
systems = vcat(system_solving, sian, other, random)

println("Running the following systems: ", map(first, systems))

timers = []
for (name, sys) in systems
    @info "Running $name.."

    for kws in [
        (_use_divmask=false, monoms=:dense),
        (_use_divmask=true, monoms=:dense),
        (_use_divmask=false, monoms=:packed),
        (_use_divmask=true, monoms=:packed),
    ]
        @info "Options:" kws
        for t in 1:2
            TimerOutputs.enable_timer!(Groebner._TIMER); reset_timer!(Groebner._TIMER);
            groebner(sys; kws...);
            if t == 2 show(Groebner._TIMER, allocations=false); println() end
        end
        push!(timers, [name, kws, copy(Groebner._TIMER)])
    end    
end

labels = ["Name", "Options", "Select Pairs", "Symb Prepr", "Reduction", "Update", "Autoreduce", "IO", "Total"]
data = []
for entry in timers
    name = entry[1]
    kws = entry[2]
    timer = entry[3]
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
=#
