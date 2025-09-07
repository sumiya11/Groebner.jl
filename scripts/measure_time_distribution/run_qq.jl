using AbstractAlgebra, Groebner, BenchmarkTools, PrettyTables, Printf, TimerOutputs

include("gizmos.jl")
include("../mqrr/mqrr.jl")

function Groebner.ratrec_nemo(a::Nemo.ZZRingElem, m::Nemo.ZZRingElem, N::Nemo.ZZRingElem, D::Nemo.ZZRingElem)
    ok, n, d = mqrr(BigInt(a), BigInt(m))
    ok, d != 0 ? Rational{BigInt}(n, d) : Rational{BigInt}(0)
end

k = AbstractAlgebra.QQ
system_solving = [
    ("chandra-9", Groebner.Examples.chandran(9, k=k)),
    ("chandra-10", Groebner.Examples.chandran(10, k=k)),
    ("chandra-11", Groebner.Examples.chandran(11, k=k)),
    # ("chandra-12", Groebner.Examples.chandran(12, k=k)),
    # ("chandra-13", Groebner.Examples.chandran(13, k=k)),
    ("cyclic-7", Groebner.Examples.cyclicn(7, k=k)),
    ("cyclic-8", Groebner.Examples.cyclicn(8, k=k)),
    # ("cyclic-9", Groebner.Examples.cyclicn(9, k=k)),
    ("eco-11", Groebner.Examples.econ(11, k=k)),
    ("eco-12", Groebner.Examples.econ(12, k=k)),
    ("eco-13", Groebner.Examples.econ(13, k=k)),
    # ("eco-14", Groebner.Examples.econ(14, k=k)),
    ("noon-7", Groebner.Examples.noonn(7, k=k)),
    ("noon-8", Groebner.Examples.noonn(8, k=k)),
    ("noon-9", Groebner.Examples.noonn(9, k=k)),
    # ("noon-10", Groebner.Examples.noonn(10, k=k)),
    ("henrion-6", Groebner.Examples.henrion6(k=k)),
    # ("henrion-7", Groebner.Examples.henrion7(k=k)),
    # ("henrion-8", Groebner.Examples.henrion8(k=k)),
    ("katsura-9", Groebner.Examples.katsuran(9, k=k)),
    ("katsura-10", Groebner.Examples.katsuran(10, k=k)),
    ("katsura-11", Groebner.Examples.katsuran(11, k=k)),
    # ("katsura-12", Groebner.Examples.katsuran(12, k=k)),
    ("reimer-6", Groebner.Examples.reimern(6, k=k)),
    ("reimer-7", Groebner.Examples.reimern(7, k=k)),
    ("reimer-8", Groebner.Examples.reimern(8, k=k)),
    ("hexapod", Groebner.Examples.hexapod(k=k)),
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

timers = []
for (name, sys) in systems
    @info "Running $name.."
    for t in 1:2
        TimerOutputs.enable_timer!(Groebner._TIMER); 
        reset_timer!(Groebner._TIMER); 
        groebner(sys); 
        if t == 2 show(Groebner._TIMER, allocations=false); println() end
    end
    push!(timers, [name, copy(Groebner._TIMER)])
end

labels = ["Name", "Guess Prime", "Learn", "Apply", "CRT", "RatRec", "Check", "IO", "Total"]
data, data2 = [], []
for entry in timers
    name = entry[1]
    timer = entry[2]
    time_guess_prime = TimerOutputs.time(timer["_groebner_guess_lucky_prime"])
    time_learn = TimerOutputs.time(timer["f4_learn!"])
    time_apply = TimerOutputs.time(timer["f4_apply!"])
    time_crt = TimerOutputs.time(timer["crt_vec_full!"]) + 
               TimerOutputs.time(timer["crt_vec_partial!"])
    time_ratrec = TimerOutputs.time(timer["ratrec_vec_full!"]) + 
                  TimerOutputs.time(timer["ratrec_vec_partial!"])
    time_io = TimerOutputs.time(timer["io_convert_polynomials_to_ir"]) +
              TimerOutputs.time(timer["io_convert_ir_to_polynomials"]) +
              TimerOutputs.time(timer["ir_convert_internal_to_ir"]) +
                TimerOutputs.time(timer["ir_convert_ir_to_internal"])
    time_check = TimerOutputs.time(timer["modular_lift_check!"])
    n_primes = TimerOutputs.ncalls(timer["f4_apply!"])
    time_total = TimerOutputs.tottime(timer)
    push!(data, [name, time_guess_prime, time_learn, time_apply, time_crt, time_ratrec, time_check, time_io, time_total])
    push!(data2, [name, 4*n_primes])
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
        (v,i,j) -> v isa Number ? (j != size(matrix, 2) ? "$(round(Int, 100*v / matrix[i,end]))%" : @sprintf("%.2f", v / 1e9)) : v,
    ],
    # row_group_labels = [1 => "System solving", length(system_solving)+1 => "SIAN", length(system_solving)+length(sian)+1 => "Other"],
    highlighters  = [hl5, hl6],
    table_format = TextTableFormat(borders = text_table_borders__ascii_rounded),
    summary_row_labels = ["Total"],
    summary_rows = [(data, i) -> i > 1 ? @sprintf("%d%%", round(Int, 100*sum(data[:, i]) / sum(data[:, end]))) : ""],
    fit_table_in_display_vertically = false,
)

matrix2 = permutedims(reduce(hcat, data2))

pretty_table(
    matrix2, 
    column_labels=[labels[1], "Primes"],
    title="Other statistics",
    table_format = TextTableFormat(borders = text_table_borders__ascii_rounded),
    summary_row_labels = ["Total"],
    summary_rows = [(data, i) -> i > 1 ? sum(data[:, i]) : ""],
    fit_table_in_display_vertically = false,
)

#=
Default output:

                                    Timings in seconds
.-------.------------.-------------.-------.-------.-----.--------.-------.-----.--------.
|       |       Name | Guess Prime | Learn | Apply | CRT | RatRec | Check |  IO |  Total |
:-------+------------+-------------+-------+-------+-----+--------+-------+-----+--------:
|       |  chandra-9 |          7% |    7% |   24% | 26% |    28% |    7% |  0% |   0.55 |
|       | chandra-10 |          4% |    5% |   55% | 18% |    14% |    4% |  0% |   3.96 |
|       | chandra-11 |          9% |   19% |   15% | 29% |    22% |    6% |  0% |   9.23 |
|       | chandra-12 |          9% |   14% |   18% | 31% |    23% |    5% |  0% |  45.31 |
|       | chandra-13 |         11% |   21% |   20% | 26% |    18% |    4% |  0% | 211.95 |
|       |   cyclic-7 |         19% |   23% |   23% |  9% |    19% |    7% |  1% |   0.54 |
|       |   cyclic-8 |         15% |   26% |   47% |  3% |     7% |    2% |  0% |  11.96 |
|       |     eco-11 |         27% |   36% |   17% |  5% |    10% |    5% |  0% |   2.54 |
|       |     eco-12 |         19% |   44% |   24% |  3% |     6% |    3% |  0% |  15.21 |
|       |     eco-13 |         14% |   55% |   18% |  4% |     7% |    1% |  1% | 128.80 |
|       |     noon-7 |          9% |    5% |    1% |  1% |     3% |    3% | 79% |   2.64 |
|       |     noon-8 |         39% |   23% |    4% |  4% |    10% |   19% |  1% |   4.75 |
|       |     noon-9 |         50% |   23% |    4% |  3% |     9% |   10% |  1% |  41.79 |
|       |  henrion-6 |         11% |   34% |   14% | 12% |    23% |    6% |  1% |   0.42 |
|       |  katsura-9 |         18% |   28% |   16% | 10% |    23% |    5% |  0% |   1.90 |
|       | katsura-10 |          8% |   41% |   21% |  9% |    14% |    6% |  0% |  15.93 |
|       | katsura-11 |          8% |   35% |   37% |  8% |    10% |    2% |  0% | 114.65 |
|       |   reimer-6 |         27% |   24% |   14% |  9% |    15% |    9% |  1% |   0.20 |
|       |   reimer-7 |         25% |   29% |   17% |  9% |    14% |    5% |  0% |   3.60 |
|       |   reimer-8 |         20% |   33% |   25% |  9% |     9% |    4% |  0% | 128.94 |
|       |    hexapod |          0% |    0% |   10% | 57% |    33% |    0% |  0% |   3.68 |
|       |        ipp |          0% |    0% |    4% | 75% |    21% |    0% |  0% |  24.53 |
|       |       hiv2 |         20% |   12% |   68% |  0% |     0% |    0% |  0% |  12.13 |
|       |      alea6 |          1% |    1% |   23% | 39% |    35% |    1% |  0% |  18.29 |
:-------+------------+-------------+-------+-------+-----+--------+-------+-----+--------:
| Total |            |         14% |   30% |   23% | 16% |    13% |    4% |  1% |   100% |
'-------'------------'-------------'-------'-------'-----'--------'-------'-----'--------'

       Other statistics
.-------.------------.--------.
|       |       Name | Primes |
:-------+------------+--------:
|       |  chandra-9 |     72 |
|       | chandra-10 |     88 |
|       | chandra-11 |    100 |
|       | chandra-12 |    136 |
|       | chandra-13 |    152 |
|       |   cyclic-7 |     20 |
|       |   cyclic-8 |     48 |
|       |     eco-11 |     12 |
|       |     eco-12 |     16 |
|       |     eco-13 |     24 |
|       |     noon-7 |      4 |
|       |     noon-8 |      4 |
|       |     noon-9 |      4 |
|       |  henrion-6 |     24 |
|       |  katsura-9 |     28 |
|       | katsura-10 |     48 |
|       | katsura-11 |     72 |
|       |   reimer-6 |     12 |
|       |   reimer-7 |     28 |
|       |   reimer-8 |     72 |
|       |    hexapod |   1004 |
|       |        ipp |   2168 |
|       |       hiv2 |    152 |
|       |      alea6 |    380 |
:-------+------------+--------:
| Total |            |   4668 |
'-------'------------'--------'

MQRR output:
                                   Timings in seconds
.-------.------------.-------------.-------.-------.-----.--------.-------.----.--------.
|       |       Name | Guess Prime | Learn | Apply | CRT | RatRec | Check | IO |  Total |
:-------+------------+-------------+-------+-------+-----+--------+-------+----+--------:
|       |  chandra-9 |          0% |    0% |    1% |  1% |    89% |    8% | 1% |  12.65 |
|       | chandra-10 |          0% |    2% |    1% |  2% |    94% |    0% | 0% |  42.84 |
|       | chandra-11 |          1% |    1% |    1% |  2% |    95% |    0% | 0% | 152.63 |
|       |   cyclic-7 |          1% |    2% |    2% |  1% |    93% |    1% | 0% |   7.29 |
|       |   cyclic-8 |          3% |    3% |    9% |  0% |    84% |    0% | 0% |  94.75 |
|       |     eco-11 |          3% |   18% |    3% |  1% |    73% |    2% | 0% |  15.30 |
|       |     eco-12 |          6% |   11% |    9% |  1% |    72% |    1% | 0% |  62.59 |
|       |     eco-13 |          2% |    7% |    3% |  0% |    87% |    1% | 0% | 993.22 |
|       |     noon-7 |          7% |    4% |    1% |  1% |    83% |    4% | 0% |   3.03 |
|       |     noon-8 |          9% |    9% |    1% |  1% |    75% |    4% | 0% |  20.49 |
|       |     noon-9 |         17% |   11% |    2% |  1% |    62% |    6% | 1% | 115.05 |
|       |  henrion-6 |          1% |    1% |    1% |  1% |    94% |    1% | 0% |   4.71 |
|       |  katsura-9 |          1% |    2% |    1% |  1% |    95% |    1% | 0% |  29.11 |
|       | katsura-10 |          1% |    2% |    3% |  1% |    92% |    1% | 0% | 171.66 |
|       | katsura-11 |          1% |    5% |    5% |  1% |    88% |    0% | 0% | 893.24 |
|       |   reimer-6 |          2% |    2% |    1% |  1% |    92% |    2% | 0% |   2.31 |
|       |   reimer-7 |          3% |    7% |    2% |  1% |    85% |    1% | 0% |  29.07 |
|       |   reimer-8 |          3% |    5% |    4% |  1% |    85% |    1% | 0% | 823.20 |
|       |    hexapod |          0% |    0% |    0% |  2% |    97% |    0% | 0% |  86.82 |
|       |        ipp |          0% |    0% |    0% |  5% |    94% |    0% | 0% | 295.62 |
|       |       hiv2 |         27% |   19% |   53% |  0% |     1% |    0% | 0% |   9.18 |
|       |      alea6 |          0% |    0% |    1% |  1% |    97% |    0% | 0% | 491.19 |
:-------+------------+-------------+-------+-------+-----+--------+-------+----+--------:
| Total |            |          2% |    4% |    3% |  1% |    88% |    1% | 0% |   100% |
'-------'------------'-------------'-------'-------'-----'--------'-------'----'--------'
       Other statistics
.-------.------------.--------.
|       |       Name | Primes |
:-------+------------+--------:
|       |  chandra-9 |     72 |
|       | chandra-10 |     88 |
|       | chandra-11 |    100 |
|       |   cyclic-7 |     20 |
|       |   cyclic-8 |     48 |
|       |     eco-11 |     12 |
|       |     eco-12 |     16 |
|       |     eco-13 |     24 |
|       |     noon-7 |      4 |
|       |     noon-8 |      4 |
|       |     noon-9 |      4 |
|       |  henrion-6 |     24 |
|       |  katsura-9 |     28 |
|       | katsura-10 |     48 |
|       | katsura-11 |     72 |
|       |   reimer-6 |     12 |
|       |   reimer-7 |     28 |
|       |   reimer-8 |     72 |
|       |    hexapod |   1004 |
|       |        ipp |   2168 |
|       |       hiv2 |     88 |
|       |      alea6 |    380 |
:-------+------------+--------:
| Total |            |   4316 |
'-------'------------'--------'


ONLY ONE SYSTEM: hiv2
=#
