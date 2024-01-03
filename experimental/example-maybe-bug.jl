using AbstractAlgebra, BenchmarkTools, Base.Threads, Profile, Primes # , Groebner

macro my_profview(ex)
    :((VSCodeServer.Profile).clear();
    VSCodeServer.Profile.init(n=10^8, delay=0.0001);
    VSCodeServer.Profile.start_timer();
    $ex;
    VSCodeServer.Profile.stop_timer();
    VSCodeServer.view_profile(;))
end

@info "" nthreads()
@show ENV["JULIA_NUM_THREADS"]
Groebner.logging_enabled() = true
Groebner.invariants_enabled() = false
Groebner.performance_counters_enabled() = false

R, (x1, x2, x3) = polynomial_ring(QQ, ["x1", "x2", "x3"], ordering=:degrevlex)
k = [
    x1^2 * x2^2 + BigInt(2)^1000 * x3^2,
    x1^3 * x2 * x3^3 - BigInt(2)^1000 - 1 * x3,
    x1 * x2 * x3 + x1 * x2 + x2 * x3 + 11
]

k = Groebner.noonn(2, ground=GF(2^31 - 1), ordering=:degrevlex)

@my_profview Groebner.groebner(k, loglevel=0);
Groebner.groebner(k, loglevel=-3);
@btime Groebner.groebner($k, loglevel=0);
@btime Groebner.groebner($k, loglevel=0, threaded=:yes);

@profile Groebner.groebner(k, use_flint=true);

#! format: off
R,(t1,t2,t3,a,b,c) = polynomial_ring(QQ, ["t1","t2","t3","a", "b", "c"], ordering=:degrevlex)
hexapod = [1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1065102000*a^2*t1-1566200000*a^2*t2+359610000*a^2*t3-4000000*a*b*t2-1574352000*a*b*t3+4000000*a*c*t1+273640000*a*c*t3-1065102000*b^2*t1+8152000*b^2*t2+355610000*b^2*t3-1574352000*b*c*t1-273640000*b*c*t2-791462000*c^2*t1-1566200000*c^2*t2+355610000*c^2*t3+740236705137*a^2-279943961360*a*b+47071636200*a*c+1574352000*a*t1-273640000*a*t2+126292488913*b^2+837307375312*b*c+4000000*b*t1-273640000*b*t3+612513941897*c^2+4000000*c*t2-1574352000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-624135247952*a-50784764200*b-283060057360*c-791462000*t1+8152000*t2+359610000*t3+165673, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1889130000*a^2*t1-139016000*a^2*t2+357608000*a^2*t3+550492000*a*b*t3+1500376000*a*c*t3-1889130000*b^2*t1-689508000*b^2*t2+357608000*b^2*t3+550492000*b*c*t1-1500376000*b*c*t2-388754000*c^2*t1-139016000*c^2*t2+357608000*c^2*t3+740396599024*a^2+98430171568*a*b+268273230304*a*c-550492000*a*t1-1500376000*a*t2+854420557476*b^2-2714848476*b*c-1500376000*b*t3-114024022072*c^2+550492000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2+624263610988*a-268273230304*b+98430171568*c-388754000*t1-689508000*t2+357608000*t3-63620, 4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2-3295636000*a^2*t1+6825304000*a^2*t2+1438448000*a^2*t3-16000000*a*b*t2+4096192000*a*b*t3+16000000*a*c*t1+4906624000*a*c*t3-3295636000*b^2*t1+2729112000*b^2*t2+1422448000*b^2*t3+4096192000*b*c*t1-4906624000*b*c*t2+1610988000*c^2*t1+6825304000*c^2*t2+1422448000*c^2*t3+2962666483625*a^2+722869290752*a*b+875649162944*a*c-4096192000*a*t1-4906624000*a*t2+513760438633*b^2-3361285532000*b*c+16000000*b*t1-4906624000*b*t3+2443184693353*c^2+16000000*c*t2+4096192000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2-2498705324448*a-879018458944*b+741978122752*c+1610988000*t1+2729112000*t2+1438448000*t3+440361,4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2+3295636000*a^2*t1+6824896000*a^2*t2+1430432000*a^2*t3+4094592000*a*b*t3-4906624000*a*c*t3+3295636000*b^2*t1+2730304000*b^2*t2+1430432000*b^2*t3+4094592000*b*c*t1+4906624000*b*c*t2-1610988000*c^2*t1+6824896000*c^2*t2+1430432000*c^2*t3+2961910911797*a^2+732129427968*a*b-877323997696*a*c-4094592000*a*t1+4906624000*a*t2+516620569397*b^2+3361357491776*b*c+4906624000*b*t3+2445290017525*c^2+4094592000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2+2499114213824*a+877323997696*b+732129427968*c-1610988000*t1+2730304000*t2+1430432000*t3-324875, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2+1889602000*a^2*t1-138926000*a^2*t2+359604000*a^2*t3-4000000*a*b*t2+550036000*a*b*t3+4000000*a*c*t1-1500228000*a*c*t3+1889602000*b^2*t1-688962000*b^2*t2+355604000*b^2*t3+550036000*b*c*t1+1500228000*b*c*t2+389374000*c^2*t1-138926000*c^2*t2+355604000*c^2*t3+740903906549*a^2+99175424872*a*b-265964790856*a*c-550036000*a*t1+1500228000*a*t2+854030749541*b^2+2874521168*b*c+4000000*b*t1+1500228000*b*t3-114557203083*c^2+4000000*c*t2+550036000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-623884900400*a+270522742856*b+97519648872*c+389374000*t1-688962000*t2+359604000*t3+55909, 250000*a^2*t1^2+250000*a^2*t2^2+250000*a^2*t3^2+250000*b^2*t1^2+250000*b^2*t2^2+250000*b^2*t3^2+250000*c^2*t1^2+250000*c^2*t2^2+250000*c^2*t3^2+266341000*a^2*t1-391502000*a^2*t2+89402000*a^2*t3-393620000*a*b*t3-68228000*a*c*t3+266341000*b^2*t1+2118000*b^2*t2+89402000*b^2*t3-393620000*b*c*t1+68228000*b*c*t2+198113000*c^2*t1-391502000*c^2*t2+89402000*c^2*t3+184958257568*a^2-70380830480*a*b-12199439312*a*c+393620000*a*t1+68228000*a*t2+31688927488*b^2-209385275032*b*c+68228000*b*t3+153269490056*c^2-393620000*c*t3+250000*t1^2+250000*t2^2+250000*t3^2+156251491928*a+12199439312*b-70380830480*c+198113000*t1+2118000*t2+89402000*t3+159976]
#! format: on

@time gb1 = Groebner.groebner(hexapod, loglevel=0, use_flint=false);

@profview gb2 = Groebner.groebner(hexapod, loglevel=0, use_flint=true, statistics=:timings);

@assert gb1 == gb2

k_zp = map(f -> map_coefficients(c -> GF(2^30 + 3)(c), f), k)

@time Groebner.groebner_learn(k_zp);
@time Groebner.groebner_learn(k_zp, threaded=:yes);

@time gb1 = Groebner.groebner(k, threaded=:no);
@time Groebner.groebner(k, threaded=:yes);
@time for i in 1:10
    gb2 = Groebner.groebner(k, threaded=:yes)
    @assert gb2 == gb1
end

gb1 == gb2

@time gb1 = Groebner.groebner(k);

@my_profview Groebner.groebner(k);

gb2 = Groebner.groebner([x1 + BigInt(2)^70 * x2], loglevel=-6);

@time gb2 = Groebner.groebner(k);

###

R, (x1, x2, x3) = polynomial_ring(GF(2^27 + 29), ["x1", "x2", "x3"], ordering=:degrevlex)
R, (x1, x2, x3) =
    polynomial_ring(GF(nextprime(2^31 - 1)), ["x1", "x2", "x3"], ordering=:degrevlex)

s = [x1 * x2^10 + 2, x1^11 * x3^3 + 3, x2^6 * x3 + 4 * x1^4 - 5]
# s = Groebner.katsuran(9, ordering=:degrevlex, ground=AbstractAlgebra.GF(2^31 + 11))

gb2 = Groebner.groebner(s, arithmetic=:signed, loglevel=-3);
gb3 = Groebner.groebner(s, arithmetic=:delayed, loglevel=-3);
gb4 = Groebner.groebner(s, arithmetic=:basic, loglevel=-3);

gb2 == gb3 == gb4

gb6 = Groebner.groebner(s, arithmetic=:signed, coeffstight=true, loglevel=-3);
@my_profview Groebner.groebner(s, arithmetic=:signed);
gb7 = Groebner.groebner(s, arithmetic=:delayed, coeffstight=true, loglevel=-3);
gb8 = Groebner.groebner(s, arithmetic=:basic, coeffstight=true, loglevel=-3);

@btime Groebner.groebner($s, arithmetic=:signed);
@btime Groebner.groebner($s, arithmetic=:delayed);
@btime Groebner.groebner($s, arithmetic=:basic);

@btime Groebner.groebner($s, arithmetic=:signed, coeffstight=true);
@btime Groebner.groebner($s, arithmetic=:delayed, coeffstight=true);
@btime Groebner.groebner($s, arithmetic=:basic, coeffstight=true);

gb2 == gb3 == gb4
gb2 == gb3 == gb4 == gb6 == gb7 == gb8

trace, gb1 = Groebner.groebner_learn(s);

@btime Groebner.groebner_apply!($trace, $s);

@btime Groebner.groebner_apply!($trace, $((s, s)));

gb1 == gb2 == gb3 == gb4

m = Groebner.DelayedArithmeticZp(UInt64, UInt64, Primes.prevprime(2^31 - 1))

@code_native debuginfo = :none Groebner.linalg_dense_row_mod_p!(UInt[1, 2, 3], m)

########

Groebner.logging_enabled() = false

p1, p2, p3, p4, p5, p6, p7, p8 = Primes.nextprimes(2^30 + 3, 8)
s1 = Groebner.katsuran(11, ordering=:degrevlex, ground=GF(p1));
s2 = Groebner.katsuran(11, ordering=:degrevlex, ground=GF(p2));
s3 = Groebner.katsuran(11, ordering=:degrevlex, ground=GF(p3));
s4 = Groebner.katsuran(11, ordering=:degrevlex, ground=GF(p4));
s5 = Groebner.katsuran(11, ordering=:degrevlex, ground=GF(p5));

trace, gb = Groebner.groebner_learn(s1);

flag, gb1 = Groebner.groebner_apply!(trace, s1);
@assert flag && Groebner.groebner(s1) == gb1

flag, (gb2, gb3, gb4, gb5) = Groebner.groebner_apply!(trace, (s2, s3, s4, s5));
@assert flag &&
        Groebner.groebner(s2) == gb2 &&
        Groebner.groebner(s3) == gb3 &&
        Groebner.groebner(s4) == gb4 &&
        Groebner.groebner(s5) == gb5

@btime Groebner.groebner_apply!($trace, $s1);
@btime Groebner.groebner_apply!($trace, $((s2, s3, s4, s5)));
@btime Groebner.groebner_apply!($trace, $((s2, s3, s4, s5, s2, s3, s4, s5)));
@btime Groebner.groebner_apply!(
    $trace,
    $((s2, s3, s4, s5, s2, s3, s4, s5, s2, s3, s4, s5, s2, s3, s4, s5))
);

@my_profview for _ in 1:5
    Groebner.groebner_apply!(trace, ((s2, s3, s4, s5)))
end

########

success, gb1 = Groebner.groebner_apply!(trace, s, loglevel=-3);

@time trace1, gb1 = Groebner.groebner_learn(s, loglevel=-3, arithmetic=:signed);
@time trace2, gb2 =
    Groebner.groebner_learn(s, loglevel=-3, arithmetic=:signed, coeffstight=true);
gb0 == gb1 == gb2

flag, (gb3, gb4, gb5, gb6) = Groebner.groebner_apply!(trace0, (s, s, s, s), loglevel=-3);

Groebner.groebner(s) == gb3 == gb4 == gb5 == gb6

@my_profview Groebner.groebner(s, arithmetic=:delayed);

flag, gb4 = Groebner.groebner_apply!(trace1, s, loglevel=-3, arithmetic=:signed);
flag, gb5 =
    Groebner.groebner_apply!(trace2, s, loglevel=-3, arithmetic=:signed, coeffstight=true);

gb0 == gb3 == gb4 == gb5

flag, (gb6, gb7) = Groebner.groebner_apply!(trace0, (s, s), loglevel=-3);
flag, (gb8, gb9, gb10, gb11) = Groebner.groebner_apply!(trace0, (s, s, s, s), loglevel=-3);
flag, (gb12, gb13, gb14, gb15, gb00, gb01, gb02, gb03) = Groebner.groebner_apply!(
    trace2,
    (s, s, s, s, s, s, s, s),
    arithmetic=:signed,
    coeffstight=true,
    loglevel=-3
);

gb0 == gb6 == gb7 == gb8 == gb9 == gb10 == gb11
gb0 == gb12 == gb13 == gb14 == gb15

@btime Groebner.groebner_apply!($trace0, $s);
@btime Groebner.groebner_apply!($trace0, $((s, s)));
@btime Groebner.groebner_apply!($trace0, $((s, s, s, s)));

@btime Groebner.groebner_apply!($trace1, $s, arithmetic=:signed);
@btime Groebner.groebner_apply!($trace1, $((s, s)), arithmetic=:signed);
@btime Groebner.groebner_apply!($trace1, $((s, s, s, s)), arithmetic=:signed);

@btime Groebner.groebner_apply!($trace2, $s, arithmetic=:signed, coeffstight=true);
@btime Groebner.groebner_apply!($trace2, $((s, s)), arithmetic=:signed, coeffstight=true);
@btime Groebner.groebner_apply!($trace2, $((s, s)), arithmetic=:signed, coeffstight=true);

@btime Groebner.groebner_apply!(
    $trace2,
    $((s, s, s, s, s, s, s, s)),
    arithmetic=:signed,
    coeffstight=true
);

@my_profview for _ in 1:10
    Groebner.groebner_apply!(trace2, (s, s, s, s), arithmetic=:signed, coeffstight=true)
end

########

@time Groebner.groebner(s);
@my_profview Groebner.groebner(s)

p1 = 2^31 - 1
p2 = 2^30 - 35
p3 = 2^29 - 3
p4 = 2^28 - 57
p5 = 2^27 - 39
p6 = 2^26 - 5
p7 = 2^25 - 39
for p in (p1, p2, p3, p4, p5, p6, p7)
    s = Groebner.katsuran(8, ordering=:degrevlex, ground=AbstractAlgebra.GF(p))
    @info "p < 2^$(floor(Int, log(2, p)+1))"
    @btime Groebner.groebner($s)
end

Groebner.invariants_enabled() = false
Groebner.logging_enabled() = false

BenchmarkTools.DEFAULT_PARAMETERS.samples = 3
systems = [
    ("kat7", Groebner.katsuran(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("kat8", Groebner.katsuran(8, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("kat9", Groebner.katsuran(9, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("kat10", Groebner.katsuran(10, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("cyc6", Groebner.cyclicn(6, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("cyc7", Groebner.cyclicn(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("cyc8", Groebner.cyclicn(8, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("eco10", Groebner.eco10(ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("eco11", Groebner.eco11(ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("eco12", Groebner.eco12(ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("noon6", Groebner.noonn(6, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("noon7", Groebner.noonn(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2)))
]
for (name, s) in systems
    @info "$name: randomized thr. / randomized / deterministic thr. / deterministic"
    gb1 = @btime Groebner.groebner($s, threaded=:yes)
    gb2 = @btime Groebner.groebner($s, threaded=:no)
    gb4 = @btime Groebner.groebner($s, linalg=:deterministic)
    @assert gb1 == gb2 == gb3 == gb4
end

Groebner.logging_enabled() = true
Groebner.invariants_enabled() = false

@time Groebner.groebner(s, linalg=:randomized, loglevel=0, threaded=:no);
@my_profview gb1 = Groebner.groebner(s, linalg=:randomized, threaded=:no);
@profview_allocs gb2 = Groebner.groebner(s, linalg=:randomized, threaded=:no) sample_rate =
    0.1 C = true;

gb1 == gb2

BenchmarkTools.DEFAULT_PARAMETERS.samples = 3
systems = [
    ("kat7", Groebner.katsuran(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("kat8", Groebner.katsuran(8, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("kat9", Groebner.katsuran(9, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("kat10", Groebner.katsuran(10, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("cyc6", Groebner.cyclicn(6, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("cyc7", Groebner.cyclicn(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("cyc8", Groebner.cyclicn(8, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("eco10", Groebner.eco10(ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("eco11", Groebner.eco11(ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("eco12", Groebner.eco12(ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("noon6", Groebner.noonn(6, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2))),
    ("noon7", Groebner.noonn(7, ordering=:degrevlex, ground=AbstractAlgebra.GF(p2)))
]
for (name, s) in systems
    @info "$name: randomized thr. / randomized / deterministic thr. / deterministic"
    gb1 = @btime Groebner.groebner($s, linalg=:randomized, threaded=:yes)
    gb2 = @btime Groebner.groebner($s, linalg=:randomized, threaded=:no)
    gb3 = @btime Groebner.groebner($s, linalg=:deterministic, threaded=:yes)
    gb4 = @btime Groebner.groebner($s, loglevel=0, linalg=:deterministic)
    @assert gb1 == gb2 == gb3 == gb4
end

@time trace, gb2 = Groebner.groebner_learn(s);

@my_profview suc, gb3 = Groebner.groebner_apply!(trace, s);

println()
ground = AbstractAlgebra.GF(1048583)
run_f4_ff_degrevlex_benchmarks(ground)
# for k = 7
# 188.876 ms (213205 allocations: 42.44 MiB)

#! format: off
R,(t1,t2,t3,a,b,c) = polynomial_ring(QQ, ["t1","t2","t3","a", "b", "c"], ordering=:degrevlex)
sys1 = [1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1065102000*a^2*t1-1566200000*a^2*t2+359610000*a^2*t3-4000000*a*b*t2-1574352000*a*b*t3+4000000*a*c*t1+273640000*a*c*t3-1065102000*b^2*t1+8152000*b^2*t2+355610000*b^2*t3-1574352000*b*c*t1-273640000*b*c*t2-791462000*c^2*t1-1566200000*c^2*t2+355610000*c^2*t3+740236705137*a^2-279943961360*a*b+47071636200*a*c+1574352000*a*t1-273640000*a*t2+126292488913*b^2+837307375312*b*c+4000000*b*t1-273640000*b*t3+612513941897*c^2+4000000*c*t2-1574352000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-624135247952*a-50784764200*b-283060057360*c-791462000*t1+8152000*t2+359610000*t3+165673, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1889130000*a^2*t1-139016000*a^2*t2+357608000*a^2*t3+550492000*a*b*t3+1500376000*a*c*t3-1889130000*b^2*t1-689508000*b^2*t2+357608000*b^2*t3+550492000*b*c*t1-1500376000*b*c*t2-388754000*c^2*t1-139016000*c^2*t2+357608000*c^2*t3+740396599024*a^2+98430171568*a*b+268273230304*a*c-550492000*a*t1-1500376000*a*t2+854420557476*b^2-2714848476*b*c-1500376000*b*t3-114024022072*c^2+550492000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2+624263610988*a-268273230304*b+98430171568*c-388754000*t1-689508000*t2+357608000*t3-63620, 4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2-3295636000*a^2*t1+6825304000*a^2*t2+1438448000*a^2*t3-16000000*a*b*t2+4096192000*a*b*t3+16000000*a*c*t1+4906624000*a*c*t3-3295636000*b^2*t1+2729112000*b^2*t2+1422448000*b^2*t3+4096192000*b*c*t1-4906624000*b*c*t2+1610988000*c^2*t1+6825304000*c^2*t2+1422448000*c^2*t3+2962666483625*a^2+722869290752*a*b+875649162944*a*c-4096192000*a*t1-4906624000*a*t2+513760438633*b^2-3361285532000*b*c+16000000*b*t1-4906624000*b*t3+2443184693353*c^2+16000000*c*t2+4096192000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2-2498705324448*a-879018458944*b+741978122752*c+1610988000*t1+2729112000*t2+1438448000*t3+440361,4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2+3295636000*a^2*t1+6824896000*a^2*t2+1430432000*a^2*t3+4094592000*a*b*t3-4906624000*a*c*t3+3295636000*b^2*t1+2730304000*b^2*t2+1430432000*b^2*t3+4094592000*b*c*t1+4906624000*b*c*t2-1610988000*c^2*t1+6824896000*c^2*t2+1430432000*c^2*t3+2961910911797*a^2+732129427968*a*b-877323997696*a*c-4094592000*a*t1+4906624000*a*t2+516620569397*b^2+3361357491776*b*c+4906624000*b*t3+2445290017525*c^2+4094592000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2+2499114213824*a+877323997696*b+732129427968*c-1610988000*t1+2730304000*t2+1430432000*t3-324875, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2+1889602000*a^2*t1-138926000*a^2*t2+359604000*a^2*t3-4000000*a*b*t2+550036000*a*b*t3+4000000*a*c*t1-1500228000*a*c*t3+1889602000*b^2*t1-688962000*b^2*t2+355604000*b^2*t3+550036000*b*c*t1+1500228000*b*c*t2+389374000*c^2*t1-138926000*c^2*t2+355604000*c^2*t3+740903906549*a^2+99175424872*a*b-265964790856*a*c-550036000*a*t1+1500228000*a*t2+854030749541*b^2+2874521168*b*c+4000000*b*t1+1500228000*b*t3-114557203083*c^2+4000000*c*t2+550036000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-623884900400*a+270522742856*b+97519648872*c+389374000*t1-688962000*t2+359604000*t3+55909, 250000*a^2*t1^2+250000*a^2*t2^2+250000*a^2*t3^2+250000*b^2*t1^2+250000*b^2*t2^2+250000*b^2*t3^2+250000*c^2*t1^2+250000*c^2*t2^2+250000*c^2*t3^2+266341000*a^2*t1-391502000*a^2*t2+89402000*a^2*t3-393620000*a*b*t3-68228000*a*c*t3+266341000*b^2*t1+2118000*b^2*t2+89402000*b^2*t3-393620000*b*c*t1+68228000*b*c*t2+198113000*c^2*t1-391502000*c^2*t2+89402000*c^2*t3+184958257568*a^2-70380830480*a*b-12199439312*a*c+393620000*a*t1+68228000*a*t2+31688927488*b^2-209385275032*b*c+68228000*b*t3+153269490056*c^2-393620000*c*t3+250000*t1^2+250000*t2^2+250000*t3^2+156251491928*a+12199439312*b-70380830480*c+198113000*t1+2118000*t2+89402000*t3+159976]

R,(t1,t2,t3,a,b,c) = polynomial_ring(GF(2^30+3), ["t1","t2","t3","a", "b", "c"], ordering=:degrevlex)
sys2 = [1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1065102000*a^2*t1-1566200000*a^2*t2+359610000*a^2*t3-4000000*a*b*t2-1574352000*a*b*t3+4000000*a*c*t1+273640000*a*c*t3-1065102000*b^2*t1+8152000*b^2*t2+355610000*b^2*t3-1574352000*b*c*t1-273640000*b*c*t2-791462000*c^2*t1-1566200000*c^2*t2+355610000*c^2*t3+740236705137*a^2-279943961360*a*b+47071636200*a*c+1574352000*a*t1-273640000*a*t2+126292488913*b^2+837307375312*b*c+4000000*b*t1-273640000*b*t3+612513941897*c^2+4000000*c*t2-1574352000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-624135247952*a-50784764200*b-283060057360*c-791462000*t1+8152000*t2+359610000*t3+165673, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1889130000*a^2*t1-139016000*a^2*t2+357608000*a^2*t3+550492000*a*b*t3+1500376000*a*c*t3-1889130000*b^2*t1-689508000*b^2*t2+357608000*b^2*t3+550492000*b*c*t1-1500376000*b*c*t2-388754000*c^2*t1-139016000*c^2*t2+357608000*c^2*t3+740396599024*a^2+98430171568*a*b+268273230304*a*c-550492000*a*t1-1500376000*a*t2+854420557476*b^2-2714848476*b*c-1500376000*b*t3-114024022072*c^2+550492000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2+624263610988*a-268273230304*b+98430171568*c-388754000*t1-689508000*t2+357608000*t3-63620, 4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2-3295636000*a^2*t1+6825304000*a^2*t2+1438448000*a^2*t3-16000000*a*b*t2+4096192000*a*b*t3+16000000*a*c*t1+4906624000*a*c*t3-3295636000*b^2*t1+2729112000*b^2*t2+1422448000*b^2*t3+4096192000*b*c*t1-4906624000*b*c*t2+1610988000*c^2*t1+6825304000*c^2*t2+1422448000*c^2*t3+2962666483625*a^2+722869290752*a*b+875649162944*a*c-4096192000*a*t1-4906624000*a*t2+513760438633*b^2-3361285532000*b*c+16000000*b*t1-4906624000*b*t3+2443184693353*c^2+16000000*c*t2+4096192000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2-2498705324448*a-879018458944*b+741978122752*c+1610988000*t1+2729112000*t2+1438448000*t3+440361,4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2+3295636000*a^2*t1+6824896000*a^2*t2+1430432000*a^2*t3+4094592000*a*b*t3-4906624000*a*c*t3+3295636000*b^2*t1+2730304000*b^2*t2+1430432000*b^2*t3+4094592000*b*c*t1+4906624000*b*c*t2-1610988000*c^2*t1+6824896000*c^2*t2+1430432000*c^2*t3+2961910911797*a^2+732129427968*a*b-877323997696*a*c-4094592000*a*t1+4906624000*a*t2+516620569397*b^2+3361357491776*b*c+4906624000*b*t3+2445290017525*c^2+4094592000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2+2499114213824*a+877323997696*b+732129427968*c-1610988000*t1+2730304000*t2+1430432000*t3-324875, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2+1889602000*a^2*t1-138926000*a^2*t2+359604000*a^2*t3-4000000*a*b*t2+550036000*a*b*t3+4000000*a*c*t1-1500228000*a*c*t3+1889602000*b^2*t1-688962000*b^2*t2+355604000*b^2*t3+550036000*b*c*t1+1500228000*b*c*t2+389374000*c^2*t1-138926000*c^2*t2+355604000*c^2*t3+740903906549*a^2+99175424872*a*b-265964790856*a*c-550036000*a*t1+1500228000*a*t2+854030749541*b^2+2874521168*b*c+4000000*b*t1+1500228000*b*t3-114557203083*c^2+4000000*c*t2+550036000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-623884900400*a+270522742856*b+97519648872*c+389374000*t1-688962000*t2+359604000*t3+55909, 250000*a^2*t1^2+250000*a^2*t2^2+250000*a^2*t3^2+250000*b^2*t1^2+250000*b^2*t2^2+250000*b^2*t3^2+250000*c^2*t1^2+250000*c^2*t2^2+250000*c^2*t3^2+266341000*a^2*t1-391502000*a^2*t2+89402000*a^2*t3-393620000*a*b*t3-68228000*a*c*t3+266341000*b^2*t1+2118000*b^2*t2+89402000*b^2*t3-393620000*b*c*t1+68228000*b*c*t2+198113000*c^2*t1-391502000*c^2*t2+89402000*c^2*t3+184958257568*a^2-70380830480*a*b-12199439312*a*c+393620000*a*t1+68228000*a*t2+31688927488*b^2-209385275032*b*c+68228000*b*t3+153269490056*c^2-393620000*c*t3+250000*t1^2+250000*t2^2+250000*t3^2+156251491928*a+12199439312*b-70380830480*c+198113000*t1+2118000*t2+89402000*t3+159976]
#! format: on

# Groebner.performance_counters_enabled() = true
@time gb = Groebner.groebner(sys2, loglevel=0);
@time gb = Groebner.groebner(sys2, loglevel=0, statistics=:timings);
Groebner.performance_counters_enabled() = false
@my_profview for _ in 1:100
    gb = Groebner.groebner(sys2, loglevel=0)
end

sys1 = Groebner.katsuran(9, ordering=:degrevlex, ground=GF(2^27 + 29))
sys2 = Groebner.cyclicn(8, ordering=:degrevlex, ground=GF(2^30 + 3))
sys3 = Groebner.cyclicn(7, ordering=:degrevlex, ground=GF(2^31 - 1))

@time gb1 = Groebner.groebner(sys1, arithmetic=:auto, loglevel=-3);
@time gb2 = Groebner.groebner(sys1, arithmetic=:signed, loglevel=-3);
@time gb3 = Groebner.groebner(sys1, arithmetic=:delayed, loglevel=-3);
@time gb4 = Groebner.groebner(sys1, arithmetic=:basic, loglevel=-3);

@btime Groebner.groebner($sys1, arithmetic=:auto);
@btime Groebner.groebner($sys1, arithmetic=:signed);
@btime Groebner.groebner($sys1, arithmetic=:delayed);

gb1 == gb2 == gb3 == gb4

Groebner.groebner_apply!(trace, sys3)[2] == Groebner.groebner(sys3)
Groebner.groebner_apply!(trace, (sys2, sys3))[2] ==
(Groebner.groebner(sys2), Groebner.groebner(sys3))

@btime Groebner.groebner_apply!($trace, $sys2);
@btime Groebner.groebner_apply!($trace, $(sys2, sys1));
@btime Groebner.groebner_apply!($trace, $(sys2, sys1, sys2, sys1));

@my_profview for _ in 1:1000
    Groebner.groebner_apply!(trace, (sys2, sys1, sys2, sys1))
end

###################

function widenn(a::Groebner.CompositeInt{N, UInt32}) where {N}
    Groebner.CompositeInt(UInt64.(a.data))
end

function linalg_vector_addmul_sparsedense_mod_p!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T},
    p2::A
) where {T, I, A <: NTuple{N, C}} where {N, C}
    @inbounds mul::T = coeffs[1]

    @fastmath @inbounds for j in 1:length(indices)
        idx = indices[j]
        a = row[idx] .- C.(mul) .* C.(coeffs[j])
        # a = a .+ ((a .>> 63) .& p2)
        row[idx] = a
    end
end

N = 8
x = [ntuple(i -> i, N) for _ in 1:1_000];
i = Int32.(collect(1:20:1000));
c = [Int32.(ntuple(i -> i, N)) for _ in 1:50];

@btime linalg_vector_addmul_sparsedense_mod_p!($x, $i, $c, $(ntuple(i -> i, N)))

@code_native debuginfo = :none linalg_vector_addmul_sparsedense_mod_p!(
    x,
    i,
    c,
    (ntuple(i -> i, N))
)

function reduce_dense_row_mod!(row::Vector{T}, arithmetic) where {T}
    @inbounds for j in 1:length(row)
        row[j] = Groebner.mod_p(row[j], arithmetic)
    end

    nothing
end

#=
    a1 = b1 + mul * c1
    x1 = _mul_high(a1, mod.multiplier)
    a1 - (x1 >>> mod.shift) * mod.divisor

    a2 = b2 + mul * c2
    x2 = _mul_high(a2, mod.multiplier)
    a2 - (x >>> mod.shift) * mod.divisor

    or

    _mul_high(b1, mod.multiplier) + mod.multiplier(mul * c1, mod.multiplier)

    64-bit mul, 64-bit mul, 64-bit add
                 < 32                                  < 32              < 32
             < 32   < 64                           < 32   < 64           < 32
    _mul_high(b1, mod.multiplier) + mod.multiplier(mul, mod.multiplier) * c1

    64-bit mul, 64-bit add, 128-bit mul:
                       < 64                  
                < 32  < 32 < 32   < 64
    x1 = _mul_high(b1 + mul * c1, mod.multiplier)

    x2 = _mul_high(b2 + mul * c2, mod.multiplier)

=#

function muladd_mod_p_1(b, mul, c, mod)
    a = b + mul * c
    x = Groebner._mul_high(a, mod.multiplier)
    a - (x >>> mod.shift) * mod.divisor
end

function muladd_mod_p_12(b, mul, c, p2)
    a = b .- mul .* c
    a = a .+ ((a .>> 63) .& p2)
    a
end

function muladd_mod_p_2(b, mul, c::T, mod) where {T}
    #
    precomp_mul = (widen(mul) * mod.multiplier) >>> (4 * sizeof(T))   # < 4 * < 8 => < 12
    precomp_mul_hi, precomp_mul_lo =
        (precomp_mul >>> (4 * sizeof(T))), (precomp_mul & ((1 << (4 * sizeof(T))) - 1))
    #
    a = mul * c
    x = precomp_mul_hi * c + ((precomp_mul_lo * c) >>> (4 * sizeof(T)))
    res = (a - (x >>> mod.shift) * mod.divisor) % typeof(c)
    res = b + res
    ifelse(res < mod.divisor, res, res - mod.divisor)
end

function muladd_mod_p_3(
    b::T,
    mul::T,
    c::T,
    shift,
    divisor,
    precomp_mul_hi,
    precomp_mul_lo
) where {T <: UInt64}
    a = mul * c
    x = precomp_mul_hi * c + ((precomp_mul_lo * c) >>> (4 * sizeof(T)))
    res = a - (x >>> shift) * divisor
    res = b + res
end

function f1(v, ins, cfs, m)
    @inbounds mul = v[1]
    @inbounds for i in 1:length(ins)
        idx = ins[i]
        v[idx] = muladd_mod_p_1(v[idx], mul, cfs[i], m)
    end
    v
end

function f12(v, ins, cfs, p2)
    @inbounds mul = v[1]
    @fastmath @inbounds for i in 1:length(ins)
        idx = ins[i]
        v[idx] = muladd_mod_p_12(v[idx], mul, cfs[i], p2)
    end
    v
end

function f2(v::Vector{T}, ins, cfs, m) where {T}
    @inbounds mul = v[1]
    ss = 4 * sizeof(T)
    precomp_mul = (widen(mul) * m.multiplier) >>> ss   # < 4 * < 8 => < 12
    precomp_mul_hi, precomp_mul_lo =
        (precomp_mul >>> ss) % T, (precomp_mul & ((1 << ss) - 1)) % T
    @inbounds for i in 1:length(ins)
        idx = ins[i]
        v[idx] = muladd_mod_p_3(
            v[idx],
            mul,
            cfs[i],
            m.shift,
            m.divisor,
            precomp_mul_hi,
            precomp_mul_lo
        )
    end
    v
end

n = 10000
k = 1000
@benchmark f1(v1, i2, v2, m) setup = begin
    v1 = rand(UInt, n)
    m = Base.MultiplicativeInverses.UnsignedMultiplicativeInverse(UInt64(2^30 + 3))
    v2 = rand(UInt, k)
    i2 = rand(Int32(1):Int32(n), k)
end
@benchmark f12(v1, i2, v2, (Int64(2^30 + 3))^2) setup = begin
    v1 = rand(Int, n)
    m = Base.MultiplicativeInverses.SignedMultiplicativeInverse(Int64(2^30 + 3))
    v2 = rand(Int, k)
    i2 = rand(Int32(1):Int32(n), k)
end

T = Int64
p = T(2^30 + 3)
p2 = p^2

b, mul, c = rand(1:typemax(T)) % p, rand(1:typemax(T)) % p, rand(1:typemax(T)) % p

xxx = muladd_mod_p_12(b, mul, c, p2)
xxx % p
mod(b - mul * c, p)

T = Int64
p = T(2^30 + 3)
m = Groebner.SignedArithmeticZp(Int64, Int32, p)

function reduce_dense_row_by_sparse_row_v2!(
    row::Vector{A},
    indices::Vector{I},
    coeffs::Vector{T},
    arithmetic::Groebner.SignedArithmeticZp{A, T}
) where {I, A, T}
    p2 = arithmetic.p2
    @inbounds mul = row[indices[1]] % T

    @fastmath @inbounds for j in 1:length(indices)
        idx = indices[j]
        a = row[idx] - A(mul) * A(coeffs[j])
        a = a + signbit(a) * p2
        row[idx] = a
    end

    nothing
end

@code_native debuginfo = :none Groebner.linalg_vector_addmul_sparsedense!(
    T[1],
    Int32[1],
    Int32[1],
    m
)

@code_native debuginfo = :none reduce_dense_row_by_sparse_row_v2!(
    T[1],
    Int32[1],
    Int32[1],
    m
)

1

#=
LBB0_2:                                # %L25
                                        # =>This Inner Loop Header: Depth=1
        movslq  %esi, %rsi
        movq    -8(%rdi,%rsi,8), %rdx
        movslq  (%r11,%rax,4), %rbx
        imulq   %r10, %rbx
        subq    %rbx, %rdx
        movq    %rdx, %rbx
        sarq    $63, %rbx
        andq    %rcx, %rbx
        addq    %rdx, %rbx
        movq    %rbx, -8(%rdi,%rsi,8)
        cmpq    %rax, %r9
        je      .LBB0_4
# %bb.3:                                # %L25.L25_crit_edge
                                        #   in Loop: Header=BB0_2 Depth=1
        movl    4(%r8,%rax,4), %esi
        incq    %rax
        jmp     .LBB0_2
.LBB0_4:                                # %L57
        popq    %rbx
        popq    %rbp
        .cfi_def_cfa %rsp, 8
        retq
=#

#=
.LBB0_2:                                # %L24
                                        # =>This Inner Loop Header: Depth=1
        movslq  %esi, %rsi
        movq    -8(%rdi,%rsi,8), %rdx
        movq    (%r11,%rax,8), %rbx
        imulq   %r10, %rbx
        subq    %rbx, %rdx
        movq    %rdx, %rbx
        sarq    $63, %rbx
        andq    %rcx, %rbx
        addq    %rdx, %rbx
        movq    %rbx, -8(%rdi,%rsi,8)
        cmpq    %rax, %r9
        je      .LBB0_4
# %bb.3:                                # %L24.L24_crit_edge
                                        #   in Loop: Header=BB0_2 Depth=1
        movl    4(%r8,%rax,4), %esi
        incq    %rax
        jmp     .LBB0_2
.LBB0_4:                   
=#

# for i in 1:(100)
b, mul, c = rand(T) % p, rand(T) % p, rand(T) % p
#
precomp_mul = (widen(mul) * m.multiplier) >>> (4 * sizeof(T))   # < 4 * < 8 => < 12
precomp_mul_hi, precomp_mul_lo =
    (precomp_mul >>> (4 * sizeof(T))) % T, (precomp_mul & ((1 << (4 * sizeof(T))) - 1)) % T
#   
xxx = muladd_mod_p_1(b, mul, c, m)
xxb = muladd_mod_p_2(b, mul, c, m)
xxc = muladd_mod_p_3(b, mul, c, m.shift, m.divisor, precomp_mul_hi, precomp_mul_lo)
@assert xxx == xxb == xxc
# end

@code_native debuginfo = :none f12([Int64(1)], [Int32(1)], [Int64(2)], p^2)

x = ntuple(_ -> rand(Int64), 4)
@code_native debuginfo = :none muladd_mod_p_12(x, x, x, p^2)

@code_llvm debuginfo = :none f2([UInt64(1)], [Int32(1)], [UInt64(2)], m)

@btime muladd_mod_p_1($b, $mul, $c, $m)

@btime muladd_mod_p_3(
    $b,
    $mul,
    $c,
    $(m.shift),
    $(m.divisor),
    $precomp_mul_hi,
    $precomp_mul_lo
)

@code_native debuginfo = :none muladd_mod_p_1(b, mul, c, m)

@code_native debuginfo = :none muladd_mod_p_3(
    b,
    mul,
    c,
    m.shift,
    m.divisor,
    precomp_mul_hi,
    precomp_mul_lo
)

a = Groebner.CompositeArithmeticZp(Groebner.CompositeInt((UInt(2^30 + 3), UInt(2^31 - 1))))
x = Groebner.CompositeInt((UInt(2^30 + 3)^2 - 1, UInt(2^31 - 1)^2 - 3))

@code_llvm debuginfo = :none Groebner.mod_p(x, a)
@code_native debuginfo = :none Groebner.mod_p(x, a)

@code_native debuginfo = :none Groebner.linalg_vector_addmul_sparsedense_mod_p!(
    [x],
    Int32[1],
    [x],
    a
)

x1 = prevprime(UInt16(32749))
a1 = Groebner.SpecializedArithmeticZp(x1)

@code_native debuginfo = :none reduce_dense_row_mod!([x1], a1)

x4 = Groebner.CompositeInt((Primes.nextprimes(UInt32(2^30), 4)...,))
a4 = Groebner.CompositeArithmeticZp(x4)

@code_native debuginfo = :none linalg_vector_addmul_sparsedense_mod_p!(
    [widenn(x4)],
    Int32[1],
    [x4],
    x4
)

@code_native debuginfo = :none Groebner.mod_p(x4, a4)

function uwu(x, y, z)
    x + y * z
end
function uwu2(x, y, z, a)
    Groebner.mod_p(x + y * z, a)
end

function uwu3(x, y, z, a) end

function owo(x::Tuple{T, T}, y::Tuple{T, T}, z::Tuple{T, T}) where {T}
    x .+ y .* z
end
@code_native debuginfo = :none uwu(x, x, x)
@code_native debuginfo = :none Groebner.muladd(x, x, x)
#=
julia_uwu_21438:                        # @julia_uwu_21438
        .cfi_startproc
# %bb.0:                                # %top
        pushq   %rbp
        .cfi_def_cfa_offset 16
        .cfi_offset %rbp, -16
        movq    %rsp, %rbp
        .cfi_def_cfa_register %rbp
        movq    %rdi, %rax
        vmovdqu (%rdx), %xmm0
        vmovdqu (%rcx), %xmm1
        vpsrlq  $32, %xmm1, %xmm2
        vpmuludq        %xmm0, %xmm2, %xmm2
        vpsrlq  $32, %xmm0, %xmm3
        vpmuludq        %xmm3, %xmm1, %xmm3
        vpaddq  %xmm2, %xmm3, %xmm2
        vpsllq  $32, %xmm2, %xmm2
        vpmuludq        %xmm0, %xmm1, %xmm0
        vpaddq  %xmm2, %xmm0, %xmm0
        vpaddq  (%rsi), %xmm0, %xmm0
        vmovdqu %xmm0, (%rdi)
        popq    %rbp
        .cfi_def_cfa %rsp, 8
        retq
.Lfunc_end0:
=#
@code_native debuginfo = :none uwu2(x, x, x, a)

@btime uwu($x, $x, $x)
@btime uwu2($x, $x, $x, $a)

const yy = (UInt(1), UInt(1))
@code_native debuginfo = :none uwu(yy, yy, yy)

# 2^27 + 29
# 1.124 ms (2609 allocations: 1.14 MiB)

@btime Groebner.groebner_apply!($trace, $sys);
# 2^31 - 1
# 1.294 ms (2612 allocations: 1.14 MiB)

@btime Groebner.groebner_apply!($trace, $sys);
# 2^62 + 135
# 2.420 ms (2609 allocations: 1.29 MiB)

#! format: on

Groebner.performance_counters_enabled() = true
@time gb = Groebner.groebner(
    sys,
    ordering=Groebner.DegRevLex(),
    statistics=:timings,
    modular=:learn_and_apply
);

Groebner.isgroebner(gb, ordering=Groebner.DegRevLex())
Groebner.normalform(gb, sys, ordering=Groebner.DegRevLex())

using BenchmarkTools
Groebner.invariants_enabled() = false
Groebner.logging_enabled() = false
Groebner.performance_counters_enabled() = false
k = Groebner.katsuran(10, ground=GF(2^30 + 3), ordering=:degrevlex)
@btime Groebner.groebner($k);
