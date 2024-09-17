using BenchmarkTools, AbstractAlgebra, PrettyTables, Groebner
using Base.Threads, Primes

arithm = [:basic, :signed]
coeffstight = [true]
prms = [2^25 + 35, 2^27 + 29, 2^28 + 3, 2^29 + 11, 2^30 + 3]

Groebner.logging_enabled() = false
Groebner.invariants_enabled() = false

#! format: off
function hexapod(p)
    R,(t1,t2,t3,a,b,c) = polynomial_ring(GF(p), ["t1","t2","t3","a", "b", "c"], internal_ordering=:degrevlex)
    [1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1065102000*a^2*t1-1566200000*a^2*t2+359610000*a^2*t3-4000000*a*b*t2-1574352000*a*b*t3+4000000*a*c*t1+273640000*a*c*t3-1065102000*b^2*t1+8152000*b^2*t2+355610000*b^2*t3-1574352000*b*c*t1-273640000*b*c*t2-791462000*c^2*t1-1566200000*c^2*t2+355610000*c^2*t3+740236705137*a^2-279943961360*a*b+47071636200*a*c+1574352000*a*t1-273640000*a*t2+126292488913*b^2+837307375312*b*c+4000000*b*t1-273640000*b*t3+612513941897*c^2+4000000*c*t2-1574352000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-624135247952*a-50784764200*b-283060057360*c-791462000*t1+8152000*t2+359610000*t3+165673, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1889130000*a^2*t1-139016000*a^2*t2+357608000*a^2*t3+550492000*a*b*t3+1500376000*a*c*t3-1889130000*b^2*t1-689508000*b^2*t2+357608000*b^2*t3+550492000*b*c*t1-1500376000*b*c*t2-388754000*c^2*t1-139016000*c^2*t2+357608000*c^2*t3+740396599024*a^2+98430171568*a*b+268273230304*a*c-550492000*a*t1-1500376000*a*t2+854420557476*b^2-2714848476*b*c-1500376000*b*t3-114024022072*c^2+550492000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2+624263610988*a-268273230304*b+98430171568*c-388754000*t1-689508000*t2+357608000*t3-63620, 4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2-3295636000*a^2*t1+6825304000*a^2*t2+1438448000*a^2*t3-16000000*a*b*t2+4096192000*a*b*t3+16000000*a*c*t1+4906624000*a*c*t3-3295636000*b^2*t1+2729112000*b^2*t2+1422448000*b^2*t3+4096192000*b*c*t1-4906624000*b*c*t2+1610988000*c^2*t1+6825304000*c^2*t2+1422448000*c^2*t3+2962666483625*a^2+722869290752*a*b+875649162944*a*c-4096192000*a*t1-4906624000*a*t2+513760438633*b^2-3361285532000*b*c+16000000*b*t1-4906624000*b*t3+2443184693353*c^2+16000000*c*t2+4096192000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2-2498705324448*a-879018458944*b+741978122752*c+1610988000*t1+2729112000*t2+1438448000*t3+440361,4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2+3295636000*a^2*t1+6824896000*a^2*t2+1430432000*a^2*t3+4094592000*a*b*t3-4906624000*a*c*t3+3295636000*b^2*t1+2730304000*b^2*t2+1430432000*b^2*t3+4094592000*b*c*t1+4906624000*b*c*t2-1610988000*c^2*t1+6824896000*c^2*t2+1430432000*c^2*t3+2961910911797*a^2+732129427968*a*b-877323997696*a*c-4094592000*a*t1+4906624000*a*t2+516620569397*b^2+3361357491776*b*c+4906624000*b*t3+2445290017525*c^2+4094592000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2+2499114213824*a+877323997696*b+732129427968*c-1610988000*t1+2730304000*t2+1430432000*t3-324875, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2+1889602000*a^2*t1-138926000*a^2*t2+359604000*a^2*t3-4000000*a*b*t2+550036000*a*b*t3+4000000*a*c*t1-1500228000*a*c*t3+1889602000*b^2*t1-688962000*b^2*t2+355604000*b^2*t3+550036000*b*c*t1+1500228000*b*c*t2+389374000*c^2*t1-138926000*c^2*t2+355604000*c^2*t3+740903906549*a^2+99175424872*a*b-265964790856*a*c-550036000*a*t1+1500228000*a*t2+854030749541*b^2+2874521168*b*c+4000000*b*t1+1500228000*b*t3-114557203083*c^2+4000000*c*t2+550036000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-623884900400*a+270522742856*b+97519648872*c+389374000*t1-688962000*t2+359604000*t3+55909, 250000*a^2*t1^2+250000*a^2*t2^2+250000*a^2*t3^2+250000*b^2*t1^2+250000*b^2*t2^2+250000*b^2*t3^2+250000*c^2*t1^2+250000*c^2*t2^2+250000*c^2*t3^2+266341000*a^2*t1-391502000*a^2*t2+89402000*a^2*t3-393620000*a*b*t3-68228000*a*c*t3+266341000*b^2*t1+2118000*b^2*t2+89402000*b^2*t3-393620000*b*c*t1+68228000*b*c*t2+198113000*c^2*t1-391502000*c^2*t2+89402000*c^2*t3+184958257568*a^2-70380830480*a*b-12199439312*a*c+393620000*a*t1+68228000*a*t2+31688927488*b^2-209385275032*b*c+68228000*b*t3+153269490056*c^2-393620000*c*t3+250000*t1^2+250000*t2^2+250000*t3^2+156251491928*a+12199439312*b-70380830480*c+198113000*t1+2118000*t2+89402000*t3+159976]
end
#! format: on

function benchmark_system(system, trials=20; kwargs...)
    timings = []
    GC.gc()
    j = 1
    gb = nothing
    while j < trials
        time = @elapsed gb = Groebner.groebner(system; kwargs...)
        push!(timings, time)
        if time > 100e-3
            trials = 3
        end
        j += 1
    end
    gb, minimum(timings)
end

function create_systems(p)
    [
        (
            "kat5",
            Groebner.Examples.katsuran(5, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "kat6",
            Groebner.Examples.katsuran(6, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "kat7",
            Groebner.Examples.katsuran(7, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "kat8",
            Groebner.Examples.katsuran(8, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "kat9",
            Groebner.Examples.katsuran(9, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "kat10",
            Groebner.Examples.katsuran(10, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "kat11",
            Groebner.Examples.katsuran(11, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        ("hen5", Groebner.Examples.henrion5(internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))),
        (
            "reim4",
            Groebner.Examples.reimern(4, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "reim5",
            Groebner.Examples.reimern(5, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "cyc4",
            Groebner.Examples.cyclicn(4, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "cyc5",
            Groebner.Examples.cyclicn(5, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "cyc6",
            Groebner.Examples.cyclicn(6, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "cyc7",
            Groebner.Examples.cyclicn(7, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "cyc8",
            Groebner.Examples.cyclicn(8, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        ("eco10", Groebner.Examples.eco10(internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))),
        ("eco11", Groebner.Examples.eco11(internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))),
        ("eco12", Groebner.Examples.eco12(internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))),
        ("eco13", Groebner.Examples.eco12(internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))),
        (
            "noon4",
            Groebner.Examples.noonn(4, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "noon5",
            Groebner.Examples.noonn(5, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "noon6",
            Groebner.Examples.noonn(6, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "noon7",
            Groebner.Examples.noonn(7, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        (
            "noon8",
            Groebner.Examples.noonn(8, internal_ordering=:degrevlex, k=AbstractAlgebra.GF(p))
        ),
        ("hexapod", hexapod(p))
    ]
end

table = Matrix{Any}(undef, (length(create_systems(2)), length(arithm) * length(prms)))
table_num = Matrix{Any}(undef, (length(create_systems(2)), length(arithm) * length(prms)))

for (j, p) in enumerate(prms)
    systems = create_systems(p)
    @assert length(systems) == size(table, 1)

    @info "p < 2^$(floor(Int, log(2, p)+1))"

    gbs = []
    for (k, ar) in enumerate(arithm)
        @info "arithmetic: $ar"

        push!(gbs, [])
        for (i, (name, sys)) in enumerate(systems)
            @info "$name.."
            gb1, ti1 = benchmark_system(sys, arithmetic=ar)
            table[i, k + (j - 1) * length(arithm)] = BenchmarkTools.prettytime(ti1 * 1e9)
            table_num[i, k + (j - 1) * length(arithm)] = ti1
            push!(gbs[end], gb1)
        end
    end

    @assert allequal(gbs)
end

ps = map(p -> "< 2^$(floor(Int, log(2, p)+1))", prms)
header = reduce(vcat, map(p -> map(a -> "$p, $a", arithm), ps))

hl = Highlighter((data, i, j) -> all(table_num[i, j] .<= table_num[i, :]), crayon"green bold")
highlighter = pretty_table(
    table,
    header=header,
    tf=tf_markdown,
    highlighters=(hl,),
    row_labels=map(first, create_systems(2))
)

#=
|         | < 2^26, basic | < 2^26, delayed | < 2^28, basic | < 2^28, delayed | < 2^29, basic | < 2^29, delayed | < 2^30, basic | < 2^30, delayed | < 2^31, basic | < 2^31, delayed |
|---------|---------------|-----------------|---------------|-----------------|---------------|-----------------|---------------|-----------------|---------------|-----------------|
|    kat5 |    730.800 μs |      788.500 μs |    617.500 μs |      566.200 μs |    633.500 μs |      613.100 μs |    657.000 μs |      706.000 μs |    615.600 μs |      702.200 μs |
|    kat6 |      1.784 ms |        1.750 ms |      1.791 ms |        1.733 ms |      1.938 ms |        1.898 ms |      1.735 ms |        2.033 ms |      1.809 ms |        2.547 ms |
|    kat7 |      6.885 ms |        6.635 ms |      6.746 ms |        6.425 ms |      6.551 ms |        7.412 ms |      6.076 ms |        9.662 ms |      6.546 ms |       15.591 ms |
|    kat8 |     28.246 ms |       24.333 ms |     27.613 ms |       24.947 ms |     25.912 ms |       32.133 ms |     25.023 ms |       59.271 ms |     28.239 ms |       93.528 ms |
|    kat9 |    142.830 ms |      121.198 ms |    135.497 ms |      119.091 ms |    137.884 ms |      197.900 ms |    122.112 ms |      441.779 ms |    165.224 ms |      705.643 ms |
|   kat10 |    812.048 ms |      671.743 ms |    830.687 ms |      677.113 ms |    810.229 ms |         1.103 s |    714.078 ms |         2.875 s |       1.094 s |         6.902 s |
|   kat11 |       4.932 s |         4.050 s |       5.196 s |         4.361 s |       4.401 s |         7.496 s |       4.428 s |        20.507 s |       5.678 s |        60.104 s |
|    hen5 |      1.881 ms |        1.657 ms |      1.908 ms |        1.585 ms |      1.757 ms |        1.784 ms |      1.771 ms |        1.699 ms |      1.862 ms |        2.018 ms |
|   reim4 |     15.376 ms |       14.751 ms |     16.662 ms |       13.430 ms |     13.555 ms |       12.792 ms |     13.168 ms |       17.269 ms |     16.295 ms |       15.651 ms |
|   reim5 |    781.772 ms |      842.125 ms |    775.096 ms |      763.221 ms |    701.023 ms |      858.490 ms |    675.663 ms |         1.456 s |    852.027 ms |         3.039 s |
|    cyc4 |    105.200 μs |      114.800 μs |    100.100 μs |       99.800 μs |    100.500 μs |      102.600 μs |    100.300 μs |      110.400 μs |    105.700 μs |      104.400 μs |
|    cyc5 |    711.300 μs |      821.000 μs |    658.600 μs |      646.300 μs |    644.500 μs |      667.300 μs |    642.400 μs |      765.000 μs |    740.600 μs |      696.200 μs |
|    cyc6 |      2.884 ms |        3.090 ms |      2.727 ms |        2.526 ms |      2.616 ms |        2.606 ms |      2.508 ms |        2.734 ms |      2.874 ms |        3.120 ms |
|    cyc7 |     73.650 ms |       70.406 ms |     72.386 ms |       67.035 ms |     75.787 ms |       72.931 ms |     71.618 ms |      106.726 ms |     79.233 ms |      158.859 ms |
|    cyc8 |       1.157 s |      981.642 ms |       1.131 s |         1.010 s |       1.394 s |         1.197 s |       1.093 s |         2.208 s |       1.530 s |         4.763 s |
|   eco10 |     57.877 ms |       61.051 ms |     55.883 ms |       55.968 ms |     63.984 ms |       68.353 ms |     54.592 ms |      132.491 ms |     68.345 ms |      300.625 ms |
|   eco11 |    304.619 ms |      306.374 ms |    313.741 ms |      317.257 ms |    499.326 ms |      493.978 ms |    299.208 ms |      965.387 ms |    441.400 ms |         2.528 s |
|   eco12 |       2.007 s |         1.968 s |       1.878 s |         2.424 s |       1.895 s |         3.359 s |       2.099 s |         9.267 s |       2.794 s |        27.133 s |
|   eco13 |       2.059 s |         1.852 s |       1.828 s |         2.824 s |       2.246 s |         3.242 s |       2.523 s |         9.233 s |       2.235 s |        28.807 s |
|   noon4 |    602.600 μs |      631.700 μs |    566.900 μs |      764.800 μs |    611.600 μs |      562.000 μs |    668.000 μs |      616.400 μs |    800.000 μs |      734.000 μs |
|   noon5 |      3.545 ms |        3.478 ms |      3.140 ms |        4.212 ms |      3.666 ms |        3.385 ms |      4.066 ms |        3.849 ms |      4.363 ms |        5.974 ms |
|   noon6 |     20.533 ms |       22.407 ms |     19.337 ms |       23.667 ms |     20.692 ms |       25.880 ms |     24.329 ms |       45.775 ms |     21.513 ms |       95.320 ms |
|   noon7 |    144.968 ms |      146.774 ms |    128.802 ms |      162.340 ms |    133.974 ms |      252.806 ms |    170.418 ms |      731.718 ms |    214.646 ms |         2.270 s |
|   noon8 |       1.166 s |         1.337 s |       1.093 s |         2.393 s |       1.687 s |         4.680 s |       1.173 s |        16.254 s |       1.564 s |        59.645 s |
| hexapod |      3.722 ms |        3.357 ms |      3.233 ms |        3.237 ms |      3.561 ms |        3.070 ms |      3.731 ms |        3.376 ms |      3.337 ms |        4.965 ms |
=#
