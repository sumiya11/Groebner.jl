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
    # ("noon-10", Groebner.Examples.noonn(10, k=k)),
    ("henrion-6", Groebner.Examples.henrion6(k=k)),
    ("henrion-7", Groebner.Examples.henrion7(k=k)),
#    ("henrion-8", Groebner.Examples.henrion8(k=k)),
    ("katsura-9", Groebner.Examples.katsuran(9, k=k)),
    ("katsura-10", Groebner.Examples.katsuran(10, k=k)),
    ("katsura-11", Groebner.Examples.katsuran(11, k=k)),
    # ("katsura-12", Groebner.Examples.katsuran(12, k=k)),
    ("reimer-6", Groebner.Examples.reimern(6, k=k)),
    ("reimer-7", Groebner.Examples.reimern(7, k=k)),
    # ("reimer-8", Groebner.Examples.reimern(8, k=k)),
    ("hexapod", Groebner.Examples.hexapod(k=k)),
]
sian = [
    ("cholera", Groebner.Examples.Cholera(k=k)),
    ("hiv2", Groebner.Examples.HIV2(k=k)),
    # ("goodwin", Groebner.Examples.Goodwin_with_weights(k=k)),
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
              @eval Groebner.linalg_nblocks_in_randomized(nup::Int, nlow::Int, ncols::Int) = (B = floor(Int, sqrt(nlow / 3)) + 1; #=println(B);=# B)
              groebner(sys; linalg=:randomized);
              opt = (block=:current,)
            elseif kws == 1
              @eval Groebner.linalg_nblocks_in_randomized(nup::Int, nlow::Int, ncols::Int) =begin
                R, M, N = Float64(nup), Float64(ncols), Float64(nlow);
                r = 1;
                B = max(1, floor(Int, N / sqrt((2*R*M*N - M*N + 2*r*M*N) / (2*r*M + r))))
                # println(B)
                B
              end
              groebner(sys; linalg=:randomized);
              opt = (block=:optimal_r_1)
            elseif kws == 3
              @eval Groebner.linalg_nblocks_in_randomized(nup::Int, nlow::Int, ncols::Int) =begin
                R, M, N = Float64(nup), Float64(ncols), Float64(nlow);
                r = N;
                B = max(1, floor(Int, N / sqrt((2*R*M*N - M*N + 2*r*M*N) / (2*r*M + r))))
                # println(B)
                B
              end
              groebner(sys; linalg=:randomized);
              opt = (block=:optimal_r_N)
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
.-------.------------.---------------------.--------------.------------.-----------.--------.------------.------.--------.
|       |       Name |             Options | Select Pairs | Symb Prepr | Reduction | Update | Autoreduce |   IO |  Total |
:-------+------------+---------------------+--------------+------------+-----------+--------+------------+------+--------:
|       |  chandra-9 |         optimal_r_1 |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 | 0.00 |   0.02 |
|       |  chandra-9 | (block = :current,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 | 0.00 |   0.01 |
|       |  chandra-9 |         optimal_r_N |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 | 0.00 |   0.01 |
|       | chandra-10 |         optimal_r_1 |         0.00 |       0.01 |      0.07 |   0.01 |       0.00 | 0.01 |   0.10 |
|       | chandra-10 | (block = :current,) |         0.00 |       0.01 |      0.03 |   0.02 |       0.00 | 0.00 |   0.06 |
|       | chandra-10 |         optimal_r_N |         0.00 |       0.01 |      0.02 |   0.01 |       0.01 | 0.00 |   0.05 |
|       | chandra-11 |         optimal_r_1 |         0.02 |       0.02 |      0.44 |   0.03 |       0.01 | 0.06 |   0.57 |
|       | chandra-11 | (block = :current,) |         0.01 |       0.02 |      0.14 |   0.04 |       0.01 | 0.03 |   0.25 |
|       | chandra-11 |         optimal_r_N |         0.05 |       0.02 |      0.29 |   0.03 |       0.01 | 0.01 |   0.41 |
|       | chandra-12 |         optimal_r_1 |         0.05 |       0.12 |      3.09 |   0.12 |       0.07 | 0.25 |   3.70 |
|       | chandra-12 | (block = :current,) |         0.05 |       0.33 |      0.77 |   0.12 |       0.02 | 0.11 |   1.40 |
|       | chandra-12 |         optimal_r_N |         0.10 |       0.28 |      0.74 |   0.12 |       0.02 | 0.08 |   1.34 |
|       | chandra-13 |         optimal_r_1 |         0.37 |       0.55 |     22.13 |   0.47 |       0.08 | 0.61 |  24.21 |
|       | chandra-13 | (block = :current,) |         0.27 |       0.71 |      4.44 |   0.48 |       0.09 | 0.38 |   6.37 |
|       | chandra-13 |         optimal_r_N |         0.21 |       0.72 |      4.44 |   0.49 |       0.14 | 0.56 |   6.56 |
|       |   cyclic-7 |         optimal_r_1 |         0.00 |       0.00 |      0.05 |   0.01 |       0.00 | 0.00 |   0.06 |
|       |   cyclic-7 | (block = :current,) |         0.00 |       0.00 |      0.02 |   0.01 |       0.00 | 0.00 |   0.03 |
|       |   cyclic-7 |         optimal_r_N |         0.00 |       0.00 |      0.02 |   0.01 |       0.00 | 0.00 |   0.04 |
|       |   cyclic-8 |         optimal_r_1 |         0.03 |       0.04 |      1.42 |   0.03 |       0.00 | 0.01 |   1.54 |
|       |   cyclic-8 | (block = :current,) |         0.03 |       0.23 |      0.38 |   0.03 |       0.00 | 0.00 |   0.68 |
|       |   cyclic-8 |         optimal_r_N |         0.04 |       0.03 |      0.42 |   0.03 |       0.00 | 0.01 |   0.54 |
|       |     eco-11 |         optimal_r_1 |         0.01 |       0.02 |      0.13 |   0.01 |       0.02 | 0.00 |   0.19 |
|       |     eco-11 | (block = :current,) |         0.01 |       0.02 |      0.08 |   0.01 |       0.02 | 0.00 |   0.14 |
|       |     eco-11 |         optimal_r_N |         0.01 |       0.02 |      0.08 |   0.01 |       0.02 | 0.17 |   0.31 |
|       |     eco-12 |         optimal_r_1 |         0.04 |       0.09 |      0.81 |   0.02 |       0.14 | 0.03 |   1.13 |
|       |     eco-12 | (block = :current,) |         0.03 |       0.09 |      0.63 |   0.02 |       0.14 | 0.04 |   0.95 |
|       |     eco-12 |         optimal_r_N |         0.03 |       0.26 |      0.44 |   0.02 |       0.14 | 0.01 |   0.90 |
|       |     eco-13 |         optimal_r_1 |         0.10 |       0.58 |      5.74 |   0.12 |       0.04 | 0.17 |   6.76 |
|       |     eco-13 | (block = :current,) |         0.10 |       0.59 |      2.97 |   0.16 |       0.04 | 0.17 |   4.03 |
|       |     eco-13 |         optimal_r_N |         0.17 |       0.58 |      2.66 |   0.08 |       0.06 | 0.13 |   3.67 |
|       |     noon-7 |         optimal_r_1 |         0.01 |       0.02 |      0.08 |   0.01 |       0.00 | 0.00 |   0.11 |
|       |     noon-7 | (block = :current,) |         0.01 |       0.03 |      0.03 |   0.01 |       0.00 | 0.00 |   0.08 |
|       |     noon-7 |         optimal_r_N |         0.01 |       0.02 |      0.04 |   0.01 |       0.00 | 0.00 |   0.07 |
|       |     noon-8 |         optimal_r_1 |         0.04 |       0.10 |      1.18 |   0.05 |       0.02 | 0.07 |   1.46 |
|       |     noon-8 | (block = :current,) |         0.03 |       0.10 |      0.28 |   0.05 |       0.02 | 0.04 |   0.51 |
|       |     noon-8 |         optimal_r_N |         0.03 |       0.14 |      0.40 |   0.05 |       0.02 | 0.18 |   0.82 |
|       |     noon-9 |         optimal_r_1 |         0.19 |       0.90 |     22.81 |   0.38 |       0.14 | 0.49 |  24.92 |
|       |     noon-9 | (block = :current,) |         0.14 |       0.67 |      2.86 |   0.58 |       0.14 | 0.29 |   4.69 |
|       |     noon-9 |         optimal_r_N |         0.14 |       0.95 |      5.12 |   0.38 |       0.15 | 0.33 |   7.06 |
|       |  henrion-6 |         optimal_r_1 |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 | 0.00 |   0.02 |
|       |  henrion-6 | (block = :current,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 | 0.00 |   0.02 |
|       |  henrion-6 |         optimal_r_N |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 | 0.00 |   0.02 |
|       |  henrion-7 |         optimal_r_1 |         0.09 |       0.36 |      1.20 |   0.01 |       0.01 | 0.07 |   1.74 |
|       |  henrion-7 | (block = :current,) |         0.04 |       0.21 |      0.56 |   0.01 |       0.01 | 0.29 |   1.12 |
|       |  henrion-7 |         optimal_r_N |         0.21 |       0.19 |      0.75 |   0.01 |       0.01 | 0.07 |   1.25 |
|       |  katsura-9 |         optimal_r_1 |         0.01 |       0.01 |      0.07 |   0.00 |       0.00 | 0.02 |   0.11 |
|       |  katsura-9 | (block = :current,) |         0.00 |       0.01 |      0.04 |   0.00 |       0.00 | 0.00 |   0.06 |
|       |  katsura-9 |         optimal_r_N |         0.00 |       0.01 |      0.04 |   0.00 |       0.00 | 0.00 |   0.06 |
|       | katsura-10 |         optimal_r_1 |         0.02 |       0.03 |      0.55 |   0.01 |       0.01 | 0.04 |   0.66 |
|       | katsura-10 | (block = :current,) |         0.02 |       0.03 |      0.27 |   0.01 |       0.01 | 0.07 |   0.41 |
|       | katsura-10 |         optimal_r_N |         0.03 |       0.03 |      0.25 |   0.01 |       0.01 | 0.04 |   0.37 |
|       | katsura-11 |         optimal_r_1 |         0.07 |       0.20 |      4.35 |   0.04 |       0.02 | 0.31 |   4.99 |
|       | katsura-11 | (block = :current,) |         0.07 |       0.16 |      1.81 |   0.04 |       0.02 | 0.37 |   2.48 |
|       | katsura-11 |         optimal_r_N |         0.07 |       0.33 |      1.78 |   0.04 |       0.02 | 0.14 |   2.38 |
|       |   reimer-6 |         optimal_r_1 |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 | 0.00 |   0.02 |
|       |   reimer-6 | (block = :current,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 | 0.00 |   0.02 |
|       |   reimer-6 |         optimal_r_N |         0.00 |       0.01 |      0.01 |   0.00 |       0.00 | 0.00 |   0.02 |
|       |   reimer-7 |         optimal_r_1 |         0.02 |       0.15 |      0.22 |   0.00 |       0.00 | 0.01 |   0.40 |
|       |   reimer-7 | (block = :current,) |         0.02 |       0.10 |      0.13 |   0.00 |       0.00 | 0.02 |   0.28 |
|       |   reimer-7 |         optimal_r_N |         0.02 |       0.10 |      0.14 |   0.00 |       0.00 | 0.02 |   0.28 |
|       |    hexapod |         optimal_r_1 |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 | 0.00 |   0.00 |
|       |    hexapod | (block = :current,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 | 0.00 |   0.00 |
|       |    hexapod |         optimal_r_N |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 | 0.00 |   0.00 |
|       |    cholera |         optimal_r_1 |         1.22 |       0.29 |     64.51 |   1.60 |       0.00 | 0.00 |  67.62 |
|       |    cholera | (block = :current,) |         1.15 |       0.31 |     13.57 |   1.62 |       0.00 | 0.00 |  16.65 |
|       |    cholera |         optimal_r_N |         1.21 |       0.32 |     14.68 |   1.54 |       0.00 | 0.00 |  17.75 |
|       |       hiv2 |         optimal_r_1 |         0.08 |       0.12 |      2.27 |   0.32 |       0.00 | 0.00 |   2.79 |
|       |       hiv2 | (block = :current,) |         0.06 |       0.12 |      0.33 |   0.47 |       0.00 | 0.00 |   0.98 |
|       |       hiv2 |         optimal_r_N |         0.07 |       0.12 |      0.36 |   0.31 |       0.00 | 0.00 |   0.86 |
|       |        crn |         optimal_r_1 |         0.01 |       0.01 |      0.21 |   0.01 |       0.00 | 0.00 |   0.25 |
|       |        crn | (block = :current,) |         0.01 |       0.01 |      0.05 |   0.01 |       0.00 | 0.00 |   0.09 |
|       |        crn |         optimal_r_N |         0.01 |       0.01 |      0.05 |   0.01 |       0.00 | 0.00 |   0.09 |
|       |      alea6 |         optimal_r_1 |         0.00 |       0.01 |      0.03 |   0.00 |       0.00 | 0.00 |   0.05 |
|       |      alea6 | (block = :current,) |         0.00 |       0.01 |      0.03 |   0.00 |       0.00 | 0.00 |   0.05 |
|       |      alea6 |         optimal_r_N |         0.00 |       0.17 |      0.03 |   0.00 |       0.00 | 0.00 |   0.21 |
|       |   jason210 |         optimal_r_1 |         0.08 |       0.71 |      1.99 |   0.21 |       0.18 | 0.01 |   3.19 |
|       |   jason210 | (block = :current,) |         0.25 |       0.69 |      0.70 |   0.08 |       0.19 | 0.01 |   1.92 |
|       |   jason210 |         optimal_r_N |         0.08 |       0.70 |      2.15 |   0.04 |       0.18 | 0.01 |   3.16 |
|       |   gametwo2 |         optimal_r_1 |         0.28 |       0.10 |      9.74 |   0.11 |       0.44 | 0.48 |  11.15 |
|       |   gametwo2 | (block = :current,) |         0.08 |       0.29 |      4.47 |   0.11 |       0.40 | 0.28 |   5.64 |
|       |   gametwo2 |         optimal_r_N |         0.28 |       0.10 |      4.38 |   0.11 |       0.40 | 0.32 |   5.59 |
|       |      yang1 |         optimal_r_1 |         0.65 |       5.00 |    152.43 |   3.47 |       0.28 | 0.01 | 161.83 |
|       |      yang1 | (block = :current,) |         0.35 |       4.93 |      2.05 |   3.43 |       0.28 | 0.01 |  11.06 |
|       |      yang1 |         optimal_r_N |         0.66 |       4.98 |      2.16 |   3.30 |       0.28 | 0.01 |  11.40 |
|       |   bayes148 |         optimal_r_1 |         0.58 |       4.90 |     43.15 |   3.90 |       0.35 | 0.01 |  52.90 |
|       |   bayes148 | (block = :current,) |         0.38 |       5.72 |      6.97 |   3.65 |       0.35 | 0.02 |  17.09 |
|       |   bayes148 |         optimal_r_N |         0.38 |       5.73 |      8.83 |   2.83 |       0.79 | 0.01 |  18.59 |
|       |     mayr42 |         optimal_r_1 |         0.32 |       1.55 |      9.71 |  23.06 |       0.23 | 0.01 |  34.88 |
|       |     mayr42 | (block = :current,) |         0.33 |       1.57 |      1.27 |  23.74 |       0.24 | 0.01 |  27.15 |
|       |     mayr42 |         optimal_r_N |         0.32 |       1.54 |      0.72 |  23.73 |       0.23 | 0.01 |  26.55 |
|       |   rand-4-4 |         optimal_r_1 |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 | 0.00 |   0.00 |
|       |   rand-4-4 | (block = :current,) |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 | 0.00 |   0.00 |
|       |   rand-4-4 |         optimal_r_N |         0.00 |       0.00 |      0.00 |   0.00 |       0.00 | 0.00 |   0.00 |
|       |   rand-4-5 |         optimal_r_1 |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 | 0.00 |   0.02 |
|       |   rand-4-5 | (block = :current,) |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 | 0.00 |   0.02 |
|       |   rand-4-5 |         optimal_r_N |         0.00 |       0.00 |      0.01 |   0.00 |       0.00 | 0.00 |   0.02 |
|       |   rand-5-4 |         optimal_r_1 |         0.00 |       0.01 |      0.06 |   0.00 |       0.00 | 0.01 |   0.08 |
|       |   rand-5-4 | (block = :current,) |         0.00 |       0.56 |      0.05 |   0.00 |       0.00 | 0.01 |   0.62 |
|       |   rand-5-4 |         optimal_r_N |         0.00 |       0.01 |      0.05 |   0.00 |       0.00 | 0.01 |   0.08 |
|       |   rand-5-5 |         optimal_r_1 |         0.02 |       0.13 |      0.81 |   0.01 |       0.16 | 0.07 |   1.20 |
|       |   rand-5-5 | (block = :current,) |         0.19 |       0.09 |      0.73 |   0.01 |       0.01 | 0.06 |   1.08 |
|       |   rand-5-5 |         optimal_r_N |         0.02 |       0.09 |      0.69 |   0.01 |       0.06 | 0.12 |   0.99 |
:-------+------------+---------------------+--------------+------------+-----------+--------+------------+------+--------:
| Total |            |                     |        12.16 |      51.17 |    446.76 | 101.89 |       6.78 | 7.24 | 626.04 |
'-------'------------'---------------------'--------------'------------'-----------'--------'------------'------'--------'
                                                           Count linalg ops
.-------.------------.---------------------.--------.----------------------------.---------------------.--------------.---------------.
|       |       Name |             Options |    % p |                     addmul |                 mul | rows reduced | rows reducers |
:-------+------------+---------------------+--------+----------------------------+---------------------+--------------+---------------:
|       |  chandra-9 |         optimal_r_1 | (0, 0) |          (194741, 6502484) |        (291, 50158) |     (847, 0) |          8232 |
|       |  chandra-9 | (block = :current,) | (0, 0) |          (116967, 3793045) |        (291, 31902) |     (907, 0) |          8232 |
|       |  chandra-9 |         optimal_r_N | (0, 0) |          (110349, 3586536) |        (291, 34939) |     (892, 0) |          8232 |
|       | chandra-10 |         optimal_r_1 | (0, 0) |         (859421, 41765062) |       (556, 189096) |    (1634, 0) |         22173 |
|       | chandra-10 | (block = :current,) | (0, 0) |         (468861, 21626165) |       (556, 111660) |    (1727, 0) |         22173 |
|       | chandra-10 |         optimal_r_N | (0, 0) |         (424468, 19282483) |       (556, 115343) |    (1705, 0) |         22173 |
|       | chandra-11 |         optimal_r_1 | (0, 0) |       (3844666, 276964286) |      (1078, 730466) |    (3191, 0) |         61550 |
|       | chandra-11 | (block = :current,) | (0, 0) |       (1938232, 125574185) |      (1078, 411037) |    (3337, 0) |         61550 |
|       | chandra-11 |         optimal_r_N | (0, 0) |       (1831652, 116702637) |      (1078, 433151) |    (3304, 0) |         61550 |
|       | chandra-12 |         optimal_r_1 | (0, 0) |     (17329945, 1877710261) |     (2113, 2858054) |    (6286, 0) |        173722 |
|       | chandra-12 | (block = :current,) | (0, 0) |       (7887100, 724197188) |     (2113, 1551182) |    (6508, 0) |        173722 |
|       | chandra-12 |         optimal_r_N | (0, 0) |       (7735276, 703945801) |     (2113, 1560552) |    (6460, 0) |        173722 |
|       | chandra-13 |         optimal_r_1 | (0, 0) |    (78586602, 12949011705) |    (4173, 11255375) |   (12455, 0) |        499786 |
|       | chandra-13 | (block = :current,) | (0, 0) |     (33534789, 4347498587) |     (4173, 5748975) |   (12797, 0) |        499786 |
|       | chandra-13 |         optimal_r_N | (0, 0) |     (33189632, 4318897517) |     (4173, 5874218) |   (12718, 0) |        499786 |
|       |   cyclic-7 |         optimal_r_1 | (0, 0) |         (324773, 35998100) |       (585, 127849) |    (1401, 0) |          8926 |
|       |   cyclic-7 | (block = :current,) | (0, 0) |         (158268, 18008616) |       (585, 123495) |    (1503, 0) |          8926 |
|       |   cyclic-7 |         optimal_r_N | (0, 0) |         (166748, 18260228) |       (585, 124566) |    (1475, 0) |          8926 |
|       |   cyclic-8 |         optimal_r_1 | (0, 0) |      (4174493, 1250938147) |     (1509, 1144674) |    (3460, 0) |         53051 |
|       |   cyclic-8 | (block = :current,) | (0, 0) |       (1435136, 433566324) |     (1509, 1081294) |    (3666, 0) |         53051 |
|       |   cyclic-8 |         optimal_r_N | (0, 0) |       (1503836, 443739415) |     (1509, 1091326) |    (3568, 0) |         53051 |
|       |     eco-11 |         optimal_r_1 | (0, 0) |        (1335784, 97579452) |       (477, 178127) |    (3832, 0) |         68548 |
|       |     eco-11 | (block = :current,) | (0, 0) |        (1136949, 81926626) |       (477, 145169) |    (3928, 0) |         68548 |
|       |     eco-11 |         optimal_r_N | (0, 0) |        (1062764, 77010779) |       (477, 149143) |    (3878, 0) |         68548 |
|       |     eco-12 |         optimal_r_1 | (0, 0) |       (5786253, 665191419) |       (876, 643140) |    (8458, 0) |        188457 |
|       |     eco-12 | (block = :current,) | (0, 0) |       (4867500, 551625580) |       (876, 518708) |    (8603, 0) |        188457 |
|       |     eco-12 |         optimal_r_N | (0, 0) |       (4511445, 506057932) |       (876, 517739) |    (8526, 0) |        188457 |
|       |     eco-13 |         optimal_r_1 | (0, 0) |     (22094527, 3589105676) |     (1627, 2364166) |    (4772, 0) |        544624 |
|       |     eco-13 | (block = :current,) | (0, 0) |     (18557873, 2867593017) |     (1627, 1820466) |    (4987, 0) |        544624 |
|       |     eco-13 |         optimal_r_N | (0, 0) |     (16794266, 2588162740) |     (1627, 1844415) |    (4882, 0) |        544624 |
|       |     noon-7 |         optimal_r_1 | (0, 0) |        (1927391, 46767203) |       (488, 434297) |    (1483, 0) |         78113 |
|       |     noon-7 | (block = :current,) | (0, 0) |         (577272, 13328386) |       (488, 220984) |    (1571, 0) |         78113 |
|       |     noon-7 |         optimal_r_N | (0, 0) |         (788323, 20894725) |       (488, 284713) |    (1516, 0) |         78113 |
|       |     noon-8 |         optimal_r_1 | (0, 0) |      (19456612, 799775254) |     (1330, 3702355) |    (4012, 0) |        382756 |
|       |     noon-8 | (block = :current,) | (0, 0) |       (4500473, 157654102) |     (1330, 1647909) |    (4179, 0) |        382756 |
|       |     noon-8 |         optimal_r_N | (0, 0) |       (6013266, 290823509) |     (1330, 2145279) |    (4075, 0) |        382756 |
|       |     noon-9 |         optimal_r_1 | (0, 0) |   (195906413, 15102854952) |    (3673, 31837602) |   (11044, 0) |       1833901 |
|       |     noon-9 | (block = :current,) | (0, 0) |     (34188366, 1891707215) |    (3673, 12564524) |   (11353, 0) |       1833901 |
|       |     noon-9 |         optimal_r_N | (0, 0) |     (47215289, 4182487495) |    (3673, 16667786) |   (11157, 0) |       1833901 |
|       |  henrion-6 |         optimal_r_1 | (0, 0) |           (74375, 6866507) |         (89, 36427) |     (281, 0) |          9156 |
|       |  henrion-6 | (block = :current,) | (0, 0) |           (51812, 5117342) |         (89, 35423) |     (308, 0) |          9156 |
|       |  henrion-6 |         optimal_r_N | (0, 0) |           (74375, 6866507) |         (89, 36427) |     (281, 0) |          9156 |
|       |  henrion-7 |         optimal_r_1 | (0, 0) |      (3335263, 1204960727) |      (414, 1207588) |    (1262, 0) |        120673 |
|       |  henrion-7 | (block = :current,) | (0, 0) |       (1678330, 608389170) |      (414, 1145067) |    (1345, 0) |        120673 |
|       |  henrion-7 |         optimal_r_N | (0, 0) |       (2049212, 695034716) |      (414, 1207044) |    (1267, 0) |        120673 |
|       |  katsura-9 |         optimal_r_1 | (0, 0) |         (738124, 49970460) |       (264, 109101) |     (815, 0) |         25572 |
|       |  katsura-9 | (block = :current,) | (0, 0) |         (500114, 41038011) |       (264, 102447) |     (875, 0) |         25572 |
|       |  katsura-9 |         optimal_r_N | (0, 0) |         (470460, 36274987) |       (264, 103508) |     (839, 0) |         25572 |
|       | katsura-10 |         optimal_r_1 | (0, 0) |       (3646344, 378958738) |       (528, 429487) |    (1610, 0) |         69492 |
|       | katsura-10 | (block = :current,) | (0, 0) |       (2112967, 279521918) |       (528, 397101) |    (1707, 0) |         69492 |
|       | katsura-10 |         optimal_r_N | (0, 0) |       (2055104, 254643494) |       (528, 405960) |    (1648, 0) |         69492 |
|       | katsura-11 |         optimal_r_1 | (0, 0) |     (16937970, 2833065352) |     (1040, 1673152) |    (3149, 0) |        179164 |
|       | katsura-11 | (block = :current,) | (0, 0) |      (9083763, 1961456047) |     (1040, 1536029) |    (3299, 0) |        179164 |
|       | katsura-11 |         optimal_r_N | (0, 0) |      (9038863, 1826767402) |     (1040, 1558034) |    (3208, 0) |        179164 |
|       |   reimer-6 |         optimal_r_1 | (0, 0) |           (85662, 7024602) |        (198, 58702) |     (514, 0) |          9495 |
|       |   reimer-6 | (block = :current,) | (0, 0) |           (61940, 5225166) |        (198, 51903) |     (552, 0) |          9495 |
|       |   reimer-6 |         optimal_r_N | (0, 0) |           (66080, 5157440) |        (198, 54797) |     (525, 0) |          9495 |
|       |   reimer-7 |         optimal_r_1 | (0, 0) |        (850471, 195249073) |       (422, 459285) |    (1105, 0) |         62998 |
|       |   reimer-7 | (block = :current,) | (0, 0) |        (531725, 118079783) |       (422, 386868) |    (1180, 0) |         62998 |
|       |   reimer-7 |         optimal_r_N | (0, 0) |        (565150, 120157205) |       (422, 433308) |    (1121, 0) |         62998 |
|       |    hexapod |         optimal_r_1 | (0, 0) |            (14908, 637873) |         (113, 7772) |     (290, 0) |          1581 |
|       |    hexapod | (block = :current,) | (0, 0) |            (13248, 550056) |         (113, 6884) |     (313, 0) |          1581 |
|       |    hexapod |         optimal_r_N | (0, 0) |            (13233, 549299) |         (113, 7176) |     (306, 0) |          1581 |
|       |    cholera |         optimal_r_1 | (0, 0) |   (178970995, 44901656033) |    (7412, 18617887) |   (14889, 0) |        297513 |
|       |    cholera | (block = :current,) | (0, 0) |    (70178440, 14743806559) |     (7412, 9064124) |   (15249, 0) |        297513 |
|       |    cholera |         optimal_r_N | (0, 0) |    (72500604, 15326467286) |     (7412, 9377161) |   (15137, 0) |        297513 |
|       |       hiv2 |         optimal_r_1 | (0, 0) |     (27803808, 1969531141) |     (3521, 3529156) |    (7131, 0) |        133457 |
|       |       hiv2 | (block = :current,) | (0, 0) |       (5016874, 253835403) |      (3521, 708496) |    (7446, 0) |        133457 |
|       |       hiv2 |         optimal_r_N | (0, 0) |       (5521134, 277812577) |      (3521, 813427) |    (7350, 0) |        133457 |
|       |        crn |         optimal_r_1 | (0, 0) |        (2394117, 91626599) |       (776, 144404) |    (1609, 0) |         17029 |
|       |        crn | (block = :current,) | (0, 0) |         (854896, 39712390) |        (776, 89640) |    (1700, 0) |         17029 |
|       |        crn |         optimal_r_N | (0, 0) |         (880109, 40522247) |        (776, 92058) |    (1675, 0) |         17029 |
|       |      alea6 |         optimal_r_1 | (0, 0) |         (335030, 27704263) |        (227, 91198) |     (957, 0) |         22978 |
|       |      alea6 | (block = :current,) | (0, 0) |         (307470, 28145967) |        (227, 88543) |    (1005, 0) |         22978 |
|       |      alea6 |         optimal_r_N | (0, 0) |         (292091, 25483960) |        (227, 89957) |     (967, 0) |         22978 |
|       |   jason210 |         optimal_r_1 | (0, 0) |      (46189107, 990809215) |      (897, 3941882) |    (2736, 0) |       1237652 |
|       |   jason210 | (block = :current,) | (0, 0) |       (9968139, 197061780) |       (897, 860743) |    (2950, 0) |       1237652 |
|       |   jason210 |         optimal_r_N | (0, 0) |      (46189107, 990809215) |      (897, 3941882) |    (2736, 0) |       1237652 |
|       |   gametwo2 |         optimal_r_1 | (0, 0) |      (7073369, 7926763010) |     (2587, 5203689) |    (6742, 0) |         34643 |
|       |   gametwo2 | (block = :current,) | (0, 0) |      (5078474, 5852613481) |     (2587, 4605771) |    (6932, 0) |         34643 |
|       |   gametwo2 |         optimal_r_N | (0, 0) |      (5070524, 5811243325) |     (2587, 4644469) |    (6903, 0) |         34643 |
|       |      yang1 |         optimal_r_1 | (0, 0) |  (223392048, 190667549756) |    (4695, 82096605) |   (14156, 0) |        540323 |
|       |      yang1 | (block = :current,) | (0, 0) |       (5469641, 430765304) |     (4695, 2125961) |   (14477, 0) |        540323 |
|       |      yang1 |         optimal_r_N | (0, 0) |       (6868634, 598024236) |     (4695, 2508885) |   (14371, 0) |        540323 |
|       |   bayes148 |         optimal_r_1 | (0, 0) |   (970353340, 12609344490) |    (3605, 24711974) |   (10852, 0) |       3564965 |
|       |   bayes148 | (block = :current,) | (0, 0) |      (26405535, 130113000) |      (3605, 816599) |   (11300, 0) |       3564965 |
|       |   bayes148 |         optimal_r_N | (0, 0) |     (101120608, 724358263) |     (3605, 2552067) |   (10945, 0) |       3564965 |
|       |     mayr42 |         optimal_r_1 | (0, 0) |    (544711573, 4097839296) |     (8490, 5096514) |   (25544, 0) |        773203 |
|       |     mayr42 | (block = :current,) | (0, 0) |        (7224564, 23144339) |      (8490, 116439) |   (26914, 0) |        773203 |
|       |     mayr42 |         optimal_r_N | (0, 0) |        (7473553, 23626957) |      (8490, 116660) |   (26813, 0) |        773203 |
|       |   rand-4-4 |         optimal_r_1 | (0, 0) |            (9675, 1309702) |         (59, 12688) |     (188, 0) |          1850 |
|       |   rand-4-4 | (block = :current,) | (0, 0) |           (11023, 1535907) |         (59, 12688) |     (200, 0) |          1850 |
|       |   rand-4-4 |         optimal_r_N | (0, 0) |            (9972, 1360984) |         (59, 12688) |     (191, 0) |          1850 |
|       |   rand-4-5 |         optimal_r_1 | (0, 0) |           (30043, 9407619) |        (107, 55959) |     (335, 0) |          3967 |
|       |   rand-4-5 | (block = :current,) | (0, 0) |          (34473, 10936484) |        (107, 55959) |     (354, 0) |          3967 |
|       |   rand-4-5 |         optimal_r_N | (0, 0) |           (30606, 9618767) |        (107, 55959) |     (339, 0) |          3967 |
|       |   rand-5-4 |         optimal_r_1 | (0, 0) |         (127117, 54457062) |       (205, 174504) |     (629, 0) |          9352 |
|       |   rand-5-4 | (block = :current,) | (0, 0) |         (144324, 62228052) |       (205, 174504) |     (663, 0) |          9352 |
|       |   rand-5-4 |         optimal_r_N | (0, 0) |         (129358, 55350665) |       (205, 174504) |     (640, 0) |          9352 |
|       |   rand-5-5 |         optimal_r_1 | (0, 0) |        (687040, 834416339) |      (472, 1228002) |    (1434, 0) |         26404 |
|       |   rand-5-5 | (block = :current,) | (0, 0) |        (748968, 904999332) |      (472, 1228002) |    (1494, 0) |         26404 |
|       |   rand-5-5 |         optimal_r_N | (0, 0) |        (691382, 830994350) |      (472, 1228002) |    (1456, 0) |         26404 |
:-------+------------+---------------------+--------+----------------------------+---------------------+--------------+---------------:
| Total |            |                     | (0, 0) | (3016883979, 383486662064) | (164691, 314244974) |  (487306, 0) |      33195918 |
'-------'------------'---------------------'--------'----------------------------'---------------------'--------------'---------------'
=#
