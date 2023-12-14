using BenchmarkTools, AbstractAlgebra, PrettyTables, Groebner
using Base.Threads

p = 2^30 + 3

@info "Using $(nthreads()) Julia threads"
@info "Using prime number of order 2^$(floor(log(2,p)))"

Groebner.logging_enabled() = false
Groebner.invariants_enabled() = false

#! format: off
R,(t1,t2,t3,a,b,c) = polynomial_ring(GF(2^27+29), ["t1","t2","t3","a", "b", "c"], ordering=:degrevlex)
hexapod = [1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1065102000*a^2*t1-1566200000*a^2*t2+359610000*a^2*t3-4000000*a*b*t2-1574352000*a*b*t3+4000000*a*c*t1+273640000*a*c*t3-1065102000*b^2*t1+8152000*b^2*t2+355610000*b^2*t3-1574352000*b*c*t1-273640000*b*c*t2-791462000*c^2*t1-1566200000*c^2*t2+355610000*c^2*t3+740236705137*a^2-279943961360*a*b+47071636200*a*c+1574352000*a*t1-273640000*a*t2+126292488913*b^2+837307375312*b*c+4000000*b*t1-273640000*b*t3+612513941897*c^2+4000000*c*t2-1574352000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-624135247952*a-50784764200*b-283060057360*c-791462000*t1+8152000*t2+359610000*t3+165673, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1889130000*a^2*t1-139016000*a^2*t2+357608000*a^2*t3+550492000*a*b*t3+1500376000*a*c*t3-1889130000*b^2*t1-689508000*b^2*t2+357608000*b^2*t3+550492000*b*c*t1-1500376000*b*c*t2-388754000*c^2*t1-139016000*c^2*t2+357608000*c^2*t3+740396599024*a^2+98430171568*a*b+268273230304*a*c-550492000*a*t1-1500376000*a*t2+854420557476*b^2-2714848476*b*c-1500376000*b*t3-114024022072*c^2+550492000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2+624263610988*a-268273230304*b+98430171568*c-388754000*t1-689508000*t2+357608000*t3-63620, 4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2-3295636000*a^2*t1+6825304000*a^2*t2+1438448000*a^2*t3-16000000*a*b*t2+4096192000*a*b*t3+16000000*a*c*t1+4906624000*a*c*t3-3295636000*b^2*t1+2729112000*b^2*t2+1422448000*b^2*t3+4096192000*b*c*t1-4906624000*b*c*t2+1610988000*c^2*t1+6825304000*c^2*t2+1422448000*c^2*t3+2962666483625*a^2+722869290752*a*b+875649162944*a*c-4096192000*a*t1-4906624000*a*t2+513760438633*b^2-3361285532000*b*c+16000000*b*t1-4906624000*b*t3+2443184693353*c^2+16000000*c*t2+4096192000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2-2498705324448*a-879018458944*b+741978122752*c+1610988000*t1+2729112000*t2+1438448000*t3+440361,4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2+3295636000*a^2*t1+6824896000*a^2*t2+1430432000*a^2*t3+4094592000*a*b*t3-4906624000*a*c*t3+3295636000*b^2*t1+2730304000*b^2*t2+1430432000*b^2*t3+4094592000*b*c*t1+4906624000*b*c*t2-1610988000*c^2*t1+6824896000*c^2*t2+1430432000*c^2*t3+2961910911797*a^2+732129427968*a*b-877323997696*a*c-4094592000*a*t1+4906624000*a*t2+516620569397*b^2+3361357491776*b*c+4906624000*b*t3+2445290017525*c^2+4094592000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2+2499114213824*a+877323997696*b+732129427968*c-1610988000*t1+2730304000*t2+1430432000*t3-324875, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2+1889602000*a^2*t1-138926000*a^2*t2+359604000*a^2*t3-4000000*a*b*t2+550036000*a*b*t3+4000000*a*c*t1-1500228000*a*c*t3+1889602000*b^2*t1-688962000*b^2*t2+355604000*b^2*t3+550036000*b*c*t1+1500228000*b*c*t2+389374000*c^2*t1-138926000*c^2*t2+355604000*c^2*t3+740903906549*a^2+99175424872*a*b-265964790856*a*c-550036000*a*t1+1500228000*a*t2+854030749541*b^2+2874521168*b*c+4000000*b*t1+1500228000*b*t3-114557203083*c^2+4000000*c*t2+550036000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-623884900400*a+270522742856*b+97519648872*c+389374000*t1-688962000*t2+359604000*t3+55909, 250000*a^2*t1^2+250000*a^2*t2^2+250000*a^2*t3^2+250000*b^2*t1^2+250000*b^2*t2^2+250000*b^2*t3^2+250000*c^2*t1^2+250000*c^2*t2^2+250000*c^2*t3^2+266341000*a^2*t1-391502000*a^2*t2+89402000*a^2*t3-393620000*a*b*t3-68228000*a*c*t3+266341000*b^2*t1+2118000*b^2*t2+89402000*b^2*t3-393620000*b*c*t1+68228000*b*c*t2+198113000*c^2*t1-391502000*c^2*t2+89402000*c^2*t3+184958257568*a^2-70380830480*a*b-12199439312*a*c+393620000*a*t1+68228000*a*t2+31688927488*b^2-209385275032*b*c+68228000*b*t3+153269490056*c^2-393620000*c*t3+250000*t1^2+250000*t2^2+250000*t3^2+156251491928*a+12199439312*b-70380830480*c+198113000*t1+2118000*t2+89402000*t3+159976]
#! format: on

systems = [
    ("kat5", Groebner.katsuran(5, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("kat6", Groebner.katsuran(6, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("kat7", Groebner.katsuran(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("kat8", Groebner.katsuran(8, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("kat9", Groebner.katsuran(9, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("kat10", Groebner.katsuran(10, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("hen5", Groebner.henrion5(ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("reim4", Groebner.reimern(4, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("reim5", Groebner.reimern(5, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("cyc4", Groebner.cyclicn(4, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("cyc5", Groebner.cyclicn(5, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("cyc6", Groebner.cyclicn(6, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("cyc7", Groebner.cyclicn(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("cyc8", Groebner.cyclicn(8, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("eco10", Groebner.eco10(ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("eco11", Groebner.eco11(ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("eco12", Groebner.eco12(ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("noon4", Groebner.noonn(4, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("noon5", Groebner.noonn(5, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("noon6", Groebner.noonn(6, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("noon7", Groebner.noonn(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("noon8", Groebner.noonn(8, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))),
    ("hexapod", hexapod)
]

table = Matrix{Any}(undef, (length(systems), 4))

function benchmark_system(system, trials=5; kwargs...)
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
    linalg #1 / linalg #1 threaded / linalg #2 (default) / linalg #2 threaded"""
    gb1, ti1 = benchmark_system(s; linalg=:deterministic, threaded=:no)
    gb2, ti2 = benchmark_system(s; linalg=:deterministic, threaded=:yes)
    gb3, ti3 = benchmark_system(s; linalg=:randomized, threaded=:no)
    gb4, ti4 = benchmark_system(s; linalg=:randomized, threaded=:yes)

    (ti1, ti2, ti3, ti4) =
        map(t -> BenchmarkTools.prettytime(t * 1e9), (ti1, ti2, ti3, ti4))
    table[i, :] .= (ti1, ti2, ti3, ti4)
    println("$ti1 / $ti2 / $ti3 / $ti4")

    @assert gb1 == gb2 == gb3 == gb4
end

pretty_table(
    table,
    header=["linalg #1", "linalg #1 threaded", "linalg #2 (default)", "linalg #2 threaded"],
    tf=tf_markdown,
    row_labels=map(first, systems)
)

#=
|         |  linalg #1 | linalg #1 threaded | linalg #2 (default) | linalg #2 threaded |
|---------|------------|--------------------|---------------------|--------------------|
|    kat5 |   1.435 ms |           1.325 ms |          547.100 μs |         855.700 μs |
|    kat6 |   2.456 ms |           3.146 ms |            1.709 ms |           2.179 ms |
|    kat7 |  12.759 ms |           9.484 ms |            6.647 ms |           6.648 ms |
|    kat8 |  76.593 ms |          43.593 ms |           25.981 ms |          22.833 ms |
|    kat9 | 534.610 ms |         262.078 ms |          126.060 ms |          93.090 ms |
|   kat10 |    4.553 s |            2.032 s |          709.215 ms |         450.237 ms |
|    hen5 |   2.047 ms |           2.835 ms |            1.776 ms |           2.538 ms |
|   reim4 |  13.356 ms |          18.100 ms |           14.044 ms |          18.906 ms |
|   reim5 | 719.158 ms |         707.738 ms |          637.063 ms |         639.548 ms |
|    cyc4 |  99.200 μs |         259.700 μs |          105.400 μs |         168.300 μs |
|    cyc5 | 681.000 μs |           1.110 ms |          619.700 μs |         897.500 μs |
|    cyc6 |   3.149 ms |           4.636 ms |            2.488 ms |           3.347 ms |
|    cyc7 | 152.065 ms |         100.933 ms |           74.816 ms |          61.060 ms |
|    cyc8 |    3.512 s |            2.064 s |             1.068 s |         768.843 ms |
|   eco10 | 131.248 ms |          82.306 ms |           54.571 ms |          50.196 ms |
|   eco11 | 915.857 ms |         511.648 ms |          283.292 ms |         230.956 ms |
|   eco12 |    7.624 s |            3.862 s |             1.805 s |            1.372 s |
|   noon4 | 540.900 μs |           1.001 ms |          565.700 μs |         861.800 μs |
|   noon5 |   3.320 ms |           4.484 ms |            3.215 ms |           3.957 ms |
|   noon6 |  19.383 ms |          19.045 ms |           18.735 ms |          18.887 ms |
|   noon7 | 129.403 ms |         112.976 ms |          123.975 ms |         111.059 ms |
|   noon8 |    1.075 s |         855.434 ms |             1.031 s |         873.569 ms |
| hexapod |   4.296 ms |           4.843 ms |            3.444 ms |           4.146 ms |
=#
