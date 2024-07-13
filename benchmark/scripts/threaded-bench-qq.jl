using BenchmarkTools, AbstractAlgebra, PrettyTables
using Base.Threads

if !isdefined(Main, :Groebner)
    using Groebner
end

@info "Using $(nthreads()) Julia threads"

Groebner.logging_enabled() = false
Groebner.invariants_enabled() = false

#! format: off
R,(t1,t2,t3,a,b,c) = polynomial_ring(QQ, ["t1","t2","t3","a", "b", "c"], internal_ordering=:degrevlex)
hexapod = [1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1065102000*a^2*t1-1566200000*a^2*t2+359610000*a^2*t3-4000000*a*b*t2-1574352000*a*b*t3+4000000*a*c*t1+273640000*a*c*t3-1065102000*b^2*t1+8152000*b^2*t2+355610000*b^2*t3-1574352000*b*c*t1-273640000*b*c*t2-791462000*c^2*t1-1566200000*c^2*t2+355610000*c^2*t3+740236705137*a^2-279943961360*a*b+47071636200*a*c+1574352000*a*t1-273640000*a*t2+126292488913*b^2+837307375312*b*c+4000000*b*t1-273640000*b*t3+612513941897*c^2+4000000*c*t2-1574352000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-624135247952*a-50784764200*b-283060057360*c-791462000*t1+8152000*t2+359610000*t3+165673, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1889130000*a^2*t1-139016000*a^2*t2+357608000*a^2*t3+550492000*a*b*t3+1500376000*a*c*t3-1889130000*b^2*t1-689508000*b^2*t2+357608000*b^2*t3+550492000*b*c*t1-1500376000*b*c*t2-388754000*c^2*t1-139016000*c^2*t2+357608000*c^2*t3+740396599024*a^2+98430171568*a*b+268273230304*a*c-550492000*a*t1-1500376000*a*t2+854420557476*b^2-2714848476*b*c-1500376000*b*t3-114024022072*c^2+550492000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2+624263610988*a-268273230304*b+98430171568*c-388754000*t1-689508000*t2+357608000*t3-63620, 4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2-3295636000*a^2*t1+6825304000*a^2*t2+1438448000*a^2*t3-16000000*a*b*t2+4096192000*a*b*t3+16000000*a*c*t1+4906624000*a*c*t3-3295636000*b^2*t1+2729112000*b^2*t2+1422448000*b^2*t3+4096192000*b*c*t1-4906624000*b*c*t2+1610988000*c^2*t1+6825304000*c^2*t2+1422448000*c^2*t3+2962666483625*a^2+722869290752*a*b+875649162944*a*c-4096192000*a*t1-4906624000*a*t2+513760438633*b^2-3361285532000*b*c+16000000*b*t1-4906624000*b*t3+2443184693353*c^2+16000000*c*t2+4096192000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2-2498705324448*a-879018458944*b+741978122752*c+1610988000*t1+2729112000*t2+1438448000*t3+440361,4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2+3295636000*a^2*t1+6824896000*a^2*t2+1430432000*a^2*t3+4094592000*a*b*t3-4906624000*a*c*t3+3295636000*b^2*t1+2730304000*b^2*t2+1430432000*b^2*t3+4094592000*b*c*t1+4906624000*b*c*t2-1610988000*c^2*t1+6824896000*c^2*t2+1430432000*c^2*t3+2961910911797*a^2+732129427968*a*b-877323997696*a*c-4094592000*a*t1+4906624000*a*t2+516620569397*b^2+3361357491776*b*c+4906624000*b*t3+2445290017525*c^2+4094592000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2+2499114213824*a+877323997696*b+732129427968*c-1610988000*t1+2730304000*t2+1430432000*t3-324875, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2+1889602000*a^2*t1-138926000*a^2*t2+359604000*a^2*t3-4000000*a*b*t2+550036000*a*b*t3+4000000*a*c*t1-1500228000*a*c*t3+1889602000*b^2*t1-688962000*b^2*t2+355604000*b^2*t3+550036000*b*c*t1+1500228000*b*c*t2+389374000*c^2*t1-138926000*c^2*t2+355604000*c^2*t3+740903906549*a^2+99175424872*a*b-265964790856*a*c-550036000*a*t1+1500228000*a*t2+854030749541*b^2+2874521168*b*c+4000000*b*t1+1500228000*b*t3-114557203083*c^2+4000000*c*t2+550036000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-623884900400*a+270522742856*b+97519648872*c+389374000*t1-688962000*t2+359604000*t3+55909, 250000*a^2*t1^2+250000*a^2*t2^2+250000*a^2*t3^2+250000*b^2*t1^2+250000*b^2*t2^2+250000*b^2*t3^2+250000*c^2*t1^2+250000*c^2*t2^2+250000*c^2*t3^2+266341000*a^2*t1-391502000*a^2*t2+89402000*a^2*t3-393620000*a*b*t3-68228000*a*c*t3+266341000*b^2*t1+2118000*b^2*t2+89402000*b^2*t3-393620000*b*c*t1+68228000*b*c*t2+198113000*c^2*t1-391502000*c^2*t2+89402000*c^2*t3+184958257568*a^2-70380830480*a*b-12199439312*a*c+393620000*a*t1+68228000*a*t2+31688927488*b^2-209385275032*b*c+68228000*b*t3+153269490056*c^2-393620000*c*t3+250000*t1^2+250000*t2^2+250000*t3^2+156251491928*a+12199439312*b-70380830480*c+198113000*t1+2118000*t2+89402000*t3+159976]
#! format: on

systems = [
    # ("kat5", Groebner.Examples.katsuran(5, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("kat6", Groebner.Examples.katsuran(6, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("kat7", Groebner.Examples.katsuran(7, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    (
        "kat8",
        Groebner.Examples.katsuran(8, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)
    ),
    (
        "kat9",
        Groebner.Examples.katsuran(9, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)
    ),
    (
        "kat10",
        Groebner.Examples.katsuran(10, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)
    ),
    # ("kat11", Groebner.Examples.katsuran(11, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("kat12", Groebner.Examples.katsuran(12, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("hen5", Groebner.Examples.henrion5(internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("reim4", Groebner.Examples.reimern(4, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("reim5", Groebner.Examples.reimern(5, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("cyc4", Groebner.Examples.cyclicn(4, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("cyc5", Groebner.Examples.cyclicn(5, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    (
        "cyc6",
        Groebner.Examples.cyclicn(6, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)
    ),
    (
        "cyc7",
        Groebner.Examples.cyclicn(7, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)
    ),
    (
        "cyc8",
        Groebner.Examples.cyclicn(8, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)
    ),
    ("eco10", Groebner.Examples.eco10(internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    ("eco11", Groebner.Examples.eco11(internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("eco12", Groebner.Examples.eco12(internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("eco13", Groebner.Examples.eco13(internal_ordering=:degrevlex, k=AbstractAlgebra.QQ))
    # ("noon4", Groebner.Examples.noonn(4, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("noon5", Groebner.Examples.noonn(5, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    # ("noon6", Groebner.Examples.noonn(6, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)),
    (
        "noon7",
        Groebner.Examples.noonn(7, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)
    ),
    (
        "noon8",
        Groebner.Examples.noonn(8, internal_ordering=:degrevlex, k=AbstractAlgebra.QQ)
    ),
    ("hexapod", hexapod)
]

table = Matrix{Any}(undef, (length(systems), 2))

function benchmark_system(system, trials=3; kwargs...)
    timings = []
    GC.gc()
    gb = nothing
    for _ in 1:trials
        time = @elapsed gb = Groebner.groebner(system; kwargs...)
        push!(timings, time)
    end
    gb, minimum(timings)
end

for (i, (name, s)) in enumerate(systems)
    @info """
    $name:
    total (default) / total (threaded)"""
    # gb1, ti1 = benchmark_system(s; linalg=:deterministic, threaded=:no)
    # gb2, ti2 = benchmark_system(s; linalg=:deterministic, threaded=:yes)
    gb3, ti3 = benchmark_system(s; linalg=:randomized, threaded=:no)
    gb4, ti4 = benchmark_system(s; linalg=:randomized, threaded=:yes)

    (ti3, ti4) = map(t -> BenchmarkTools.prettytime(t * 1e9), (ti3, ti4))
    # table[i, :] .= (ti1, ti2, ti3, ti4)
    table[i, :] .= (ti3, ti4)
    println("$ti3 / $ti4")

    # @assert gb1 == gb2 == gb3 == gb4
    @assert gb3 == gb4
end

pretty_table(
    table,
    header=["total (default)", "total (threaded)"],
    tf=tf_markdown,
    row_labels=map(first, systems)
)

#=
|         |  linalg #1 | linalg #1 threaded | linalg #2 (default) | linalg #2 threaded |
|---------|------------|--------------------|---------------------|--------------------|
|    kat5 | 649.900 μs |           1.170 ms |          594.000 μs |         881.500 μs |
|    kat6 |   2.537 ms |           3.185 ms |            1.711 ms |           2.243 ms |
|    kat7 |  15.282 ms |          12.524 ms |            8.027 ms |           6.787 ms |
|    kat8 |  78.290 ms |          44.962 ms |           26.484 ms |          21.662 ms |
|    kat9 | 534.722 ms |         260.607 ms |          126.397 ms |          94.352 ms |
|   kat10 |    4.623 s |            2.205 s |          756.419 ms |         472.633 ms |
|    hen5 |   2.232 ms |           2.883 ms |            1.948 ms |           2.209 ms |
|   reim4 |  15.097 ms |          18.618 ms |           12.829 ms |          15.887 ms |
|   reim5 | 761.353 ms |         730.773 ms |          672.574 ms |         654.266 ms |
|    cyc4 | 118.600 μs |         366.000 μs |          119.800 μs |         220.800 μs |
|    cyc5 | 665.000 μs |           1.187 ms |          663.300 μs |         958.200 μs |
|    cyc6 |   2.980 ms |           5.121 ms |            2.595 ms |           3.339 ms |
|    cyc7 | 158.597 ms |         117.256 ms |           78.787 ms |          66.732 ms |
|    cyc8 |    3.498 s |            2.074 s |             1.080 s |         794.454 ms |
|   eco10 | 130.946 ms |          82.692 ms |           54.565 ms |          50.291 ms |
|   eco11 | 917.577 ms |         506.479 ms |          281.912 ms |         233.156 ms |
|   eco12 |    7.901 s |            3.888 s |             1.801 s |            1.360 s |
|   noon4 | 572.000 μs |           1.193 ms |          609.300 μs |         955.800 μs |
|   noon5 |   3.485 ms |           4.406 ms |            3.247 ms |           4.047 ms |
|   noon6 |  19.640 ms |          20.017 ms |           19.235 ms |          19.568 ms |
|   noon7 | 132.158 ms |         111.749 ms |          126.367 ms |         114.704 ms |
|   noon8 |    1.069 s |         858.967 ms |             1.066 s |         894.792 ms |
| hexapod |   4.124 ms |           4.974 ms |            3.343 ms |           4.365 ms |
=#
