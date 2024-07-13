using BenchmarkTools, AbstractAlgebra, PrettyTables # Groebner
using Base.Threads

Groebner.logging_enabled() = false
Groebner.invariants_enabled() = false

#! format: off
R,(t1,t2,t3,a,b,c) = polynomial_ring(QQ, ["t1","t2","t3","a", "b", "c"], internal_ordering=:degrevlex)
hexapod = [1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1065102000*a^2*t1-1566200000*a^2*t2+359610000*a^2*t3-4000000*a*b*t2-1574352000*a*b*t3+4000000*a*c*t1+273640000*a*c*t3-1065102000*b^2*t1+8152000*b^2*t2+355610000*b^2*t3-1574352000*b*c*t1-273640000*b*c*t2-791462000*c^2*t1-1566200000*c^2*t2+355610000*c^2*t3+740236705137*a^2-279943961360*a*b+47071636200*a*c+1574352000*a*t1-273640000*a*t2+126292488913*b^2+837307375312*b*c+4000000*b*t1-273640000*b*t3+612513941897*c^2+4000000*c*t2-1574352000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-624135247952*a-50784764200*b-283060057360*c-791462000*t1+8152000*t2+359610000*t3+165673, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1889130000*a^2*t1-139016000*a^2*t2+357608000*a^2*t3+550492000*a*b*t3+1500376000*a*c*t3-1889130000*b^2*t1-689508000*b^2*t2+357608000*b^2*t3+550492000*b*c*t1-1500376000*b*c*t2-388754000*c^2*t1-139016000*c^2*t2+357608000*c^2*t3+740396599024*a^2+98430171568*a*b+268273230304*a*c-550492000*a*t1-1500376000*a*t2+854420557476*b^2-2714848476*b*c-1500376000*b*t3-114024022072*c^2+550492000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2+624263610988*a-268273230304*b+98430171568*c-388754000*t1-689508000*t2+357608000*t3-63620, 4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2-3295636000*a^2*t1+6825304000*a^2*t2+1438448000*a^2*t3-16000000*a*b*t2+4096192000*a*b*t3+16000000*a*c*t1+4906624000*a*c*t3-3295636000*b^2*t1+2729112000*b^2*t2+1422448000*b^2*t3+4096192000*b*c*t1-4906624000*b*c*t2+1610988000*c^2*t1+6825304000*c^2*t2+1422448000*c^2*t3+2962666483625*a^2+722869290752*a*b+875649162944*a*c-4096192000*a*t1-4906624000*a*t2+513760438633*b^2-3361285532000*b*c+16000000*b*t1-4906624000*b*t3+2443184693353*c^2+16000000*c*t2+4096192000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2-2498705324448*a-879018458944*b+741978122752*c+1610988000*t1+2729112000*t2+1438448000*t3+440361,4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2+3295636000*a^2*t1+6824896000*a^2*t2+1430432000*a^2*t3+4094592000*a*b*t3-4906624000*a*c*t3+3295636000*b^2*t1+2730304000*b^2*t2+1430432000*b^2*t3+4094592000*b*c*t1+4906624000*b*c*t2-1610988000*c^2*t1+6824896000*c^2*t2+1430432000*c^2*t3+2961910911797*a^2+732129427968*a*b-877323997696*a*c-4094592000*a*t1+4906624000*a*t2+516620569397*b^2+3361357491776*b*c+4906624000*b*t3+2445290017525*c^2+4094592000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2+2499114213824*a+877323997696*b+732129427968*c-1610988000*t1+2730304000*t2+1430432000*t3-324875, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2+1889602000*a^2*t1-138926000*a^2*t2+359604000*a^2*t3-4000000*a*b*t2+550036000*a*b*t3+4000000*a*c*t1-1500228000*a*c*t3+1889602000*b^2*t1-688962000*b^2*t2+355604000*b^2*t3+550036000*b*c*t1+1500228000*b*c*t2+389374000*c^2*t1-138926000*c^2*t2+355604000*c^2*t3+740903906549*a^2+99175424872*a*b-265964790856*a*c-550036000*a*t1+1500228000*a*t2+854030749541*b^2+2874521168*b*c+4000000*b*t1+1500228000*b*t3-114557203083*c^2+4000000*c*t2+550036000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-623884900400*a+270522742856*b+97519648872*c+389374000*t1-688962000*t2+359604000*t3+55909, 250000*a^2*t1^2+250000*a^2*t2^2+250000*a^2*t3^2+250000*b^2*t1^2+250000*b^2*t2^2+250000*b^2*t3^2+250000*c^2*t1^2+250000*c^2*t2^2+250000*c^2*t3^2+266341000*a^2*t1-391502000*a^2*t2+89402000*a^2*t3-393620000*a*b*t3-68228000*a*c*t3+266341000*b^2*t1+2118000*b^2*t2+89402000*b^2*t3-393620000*b*c*t1+68228000*b*c*t2+198113000*c^2*t1-391502000*c^2*t2+89402000*c^2*t3+184958257568*a^2-70380830480*a*b-12199439312*a*c+393620000*a*t1+68228000*a*t2+31688927488*b^2-209385275032*b*c+68228000*b*t3+153269490056*c^2-393620000*c*t3+250000*t1^2+250000*t2^2+250000*t3^2+156251491928*a+12199439312*b-70380830480*c+198113000*t1+2118000*t2+89402000*t3+159976]
#! format: on

systems = [
    # ("kat5", Groebner.Examples.katsuran(5, internal_ordering=:degrevlex, k=QQ)),
    ("kat6", Groebner.Examples.katsuran(6, internal_ordering=:degrevlex, k=QQ)),
    ("kat7", Groebner.Examples.katsuran(7, internal_ordering=:degrevlex, k=QQ)),
    ("kat8", Groebner.Examples.katsuran(8, internal_ordering=:degrevlex, k=QQ)),
    ("kat9", Groebner.Examples.katsuran(9, internal_ordering=:degrevlex, k=QQ)),
    ("kat10", Groebner.Examples.katsuran(10, internal_ordering=:degrevlex, k=QQ)),
    # ("kat11", Groebner.Examples.katsuran(11, internal_ordering=:degrevlex, k=QQ)),
    # ("kat12", Groebner.Examples.katsuran(12, internal_ordering=:degrevlex, k=QQ)),
    ("hen5", Groebner.Examples.henrion5(internal_ordering=:degrevlex, k=QQ)),
    ("hen6", Groebner.Examples.henrion6(internal_ordering=:degrevlex, k=QQ)),
    ("reim4", Groebner.Examples.reimern(4, internal_ordering=:degrevlex, k=QQ)),
    # ("reim5", Groebner.Examples.reimern(5, internal_ordering=:degrevlex, k=QQ)),
    # ("cyc4", Groebner.Examples.cyclicn(4, internal_ordering=:degrevlex, k=QQ)),
    # ("cyc5", Groebner.Examples.cyclicn(5, internal_ordering=:degrevlex, k=QQ)),
    ("cyc6", Groebner.Examples.cyclicn(6, internal_ordering=:degrevlex, k=QQ)),
    ("cyc7", Groebner.Examples.cyclicn(7, internal_ordering=:degrevlex, k=QQ)),
    ("cyc8", Groebner.Examples.cyclicn(8, internal_ordering=:degrevlex, k=QQ)),
    ("eco10", Groebner.Examples.eco10(internal_ordering=:degrevlex, k=QQ)),
    ("eco11", Groebner.Examples.eco11(internal_ordering=:degrevlex, k=QQ)),
    ("eco12", Groebner.Examples.eco12(internal_ordering=:degrevlex, k=QQ)),
    # ("eco13", Groebner.Examples.eco13(internal_ordering=:degrevlex, k=QQ)),
    # ("noon4", Groebner.Examples.noonn(4, internal_ordering=:degrevlex, k=QQ)),
    ("noon5", Groebner.Examples.noonn(5, internal_ordering=:degrevlex, k=QQ)),
    ("noon6", Groebner.Examples.noonn(6, internal_ordering=:degrevlex, k=QQ)),
    ("noon7", Groebner.Examples.noonn(7, internal_ordering=:degrevlex, k=QQ)),
    ("noon8", Groebner.Examples.noonn(8, internal_ordering=:degrevlex, k=QQ)),
    ("hexapod", hexapod)
]

table = Matrix{Any}(undef, (length(systems), 3))

function benchmark_system(system, trials=20; kwargs...)
    timings = []
    GC.gc()
    j = 1
    gb = nothing
    while j < trials
        time = @elapsed gb = Groebner.groebner(system; kwargs...)
        push!(timings, time)
        if time > 100e-3
            trials = 5
        end
        if time > 10
            trials = 1
        end
        j += 1
    end
    gb, minimum(timings)
end

for (i, (name, s)) in enumerate(systems)
    @info """
    $name:
    classic multi-modular / learn & apply / learn & apply N=4"""
    gb3, ti3 = benchmark_system(s; modular=:classic_modular)
    gb4, ti4 = benchmark_system(s; modular=:learn_and_apply, batched=false)
    gb5, ti5 = benchmark_system(s; modular=:learn_and_apply, batched=true)

    (ti3, ti4, ti5) = map(t -> BenchmarkTools.prettytime(t * 1e9), (ti3, ti4, ti5))

    table[i, :] .= (ti3, ti4, ti5)
    println("$ti3 / $ti4 / $ti5")

    @assert gb3 == gb4 == gb5
end

pretty_table(
    table,
    header=["classic multi-modular", "learn & apply", "learn & apply N=4"],
    tf=tf_markdown,
    row_labels=map(first, systems)
)

#=
|         | classic multi-modular | learn & apply | learn & apply N=4 |
|---------|-----------------------|---------------|-------------------|
|    kat6 |             37.169 ms |     22.459 ms |         18.573 ms |
|    kat7 |            140.579 ms |     97.749 ms |         81.444 ms |
|    kat8 |            987.475 ms |    633.818 ms |        506.742 ms |
|    kat9 |               7.434 s |       4.626 s |           3.931 s |
|   kat10 |              76.061 s |      37.614 s |          25.093 s |
|    hen5 |            470.852 ms |    395.549 ms |        328.494 ms |
|    hen6 |               1.086 s |    655.668 ms |        522.634 ms |
|   reim4 |            151.234 ms |    106.118 ms |         81.168 ms |
|    cyc6 |             12.657 ms |      8.853 ms |          8.884 ms |
|    cyc7 |               2.255 s |       1.245 s |        884.418 ms |
|    cyc8 |              86.251 s |      47.330 s |          22.606 s |
|   eco10 |            559.307 ms |    426.300 ms |        394.343 ms |
|   eco11 |               5.185 s |       3.924 s |           2.952 s |
|   eco12 |              31.792 s |      27.085 s |          19.377 s |
|   noon5 |              9.839 ms |      8.558 ms |          8.236 ms |
|   noon6 |            114.301 ms |     65.704 ms |         65.833 ms |
|   noon7 |            653.635 ms |    392.572 ms |        396.666 ms |
|   noon8 |               5.148 s |       2.725 s |           2.816 s |
| hexapod |              15.053 s |      13.172 s |          13.193 s |
=#
