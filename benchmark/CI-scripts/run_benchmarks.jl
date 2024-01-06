# This script is not intended to be run directly.
# Instead, see runtests.jl.

using Groebner
using AbstractAlgebra
import Primes
import Nemo

suite = []

function compute_gb(system, trials=7; kws...)
    times = []
    for _ in 1:trials
        GC.gc()
        time = @elapsed groebner(system; kws...)
        push!(times, time)
    end
    minimum(times)
end

# Compute Groebner bases over integers modulo a large prime
problem = (
    problem_name="groebner, AA, GF(2^31-1), katsura 5",
    result=compute_gb(Groebner.katsuran(5, ordering=:degrevlex, ground=GF(2^31 - 1)))
)
push!(suite, problem)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), katsura 5",
        result=compute_gb(Groebner.katsuran(5, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), katsura 6",
        result=compute_gb(Groebner.katsuran(6, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), katsura 8",
        result=compute_gb(Groebner.katsuran(8, ordering=:degrevlex, ground=GF(2^31 - 1)))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), katsura 10",
        result=compute_gb(
            Groebner.katsuran(10, ordering=:degrevlex, ground=GF(2^31 - 1)),
            5
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^27+29), katsura 10",
        result=compute_gb(
            Groebner.katsuran(10, ordering=:degrevlex, ground=GF(2^27 + 29)),
            5
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^27+29), cyclic 8",
        result=compute_gb(
            Groebner.cyclicn(8, ordering=:degrevlex, ground=GF(2^27 + 29)),
            5
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, GF(2^31-1), cyclic 8",
        result=compute_gb(Groebner.cyclicn(8, ordering=:degrevlex, ground=GF(2^31 - 1)), 5)
    )
)
push!(
    suite,
    (
        problem_name="groebner, Nemo, GF(2^31-1), cyclic 8",
        result=compute_gb(
            Groebner.cyclicn(8, ordering=:degrevlex, ground=Nemo.GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner, threaded, AA, GF(2^31-1), cyclic 8",
        result=compute_gb(
            Groebner.cyclicn(8, ordering=:degrevlex, ground=GF(2^31 - 1)),
            5,
            threaded=:yes
        )
    )
)

function learn_and_apply(system)
    times = []
    trials = 7
    graph, gb = groebner_learn(system)
    for _ in 1:trials
        GC.gc()
        time2 = @elapsed groebner_apply!(graph, system)
        push!(times, time2)
    end
    minimum(times)
end

# Learn and apply
push!(
    suite,
    (
        problem_name="groebner_apply!, AA, GF(2^31-1), cyclic 7",
        result=learn_and_apply(
            Groebner.cyclicn(7, ordering=:degrevlex, ground=GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner_apply!, Nemo, GF(2^31-1), cyclic 7",
        result=learn_and_apply(
            Groebner.cyclicn(7, ordering=:degrevlex, ground=Nemo.GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner_apply!, AA, GF(2^31-1), katsura 10",
        result=learn_and_apply(
            Groebner.katsuran(10, ordering=:degrevlex, ground=GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner_apply!, AA, GF(2^27+29), katsura 10",
        result=learn_and_apply(
            Groebner.katsuran(10, ordering=:degrevlex, ground=GF(2^27 + 29))
        )
    )
)
push!(
    suite,
    (
        problem_name="groebner_apply!, Nemo, GF(2^31-1), katsura 10",
        result=learn_and_apply(
            Groebner.katsuran(10, ordering=:degrevlex, ground=Nemo.GF(2^31 - 1))
        )
    )
)

# Compute Groebner bases over the rationals
push!(
    suite,
    (
        problem_name="groebner, AA, QQ, katsura 8",
        result=compute_gb(Groebner.katsuran(8, ordering=:degrevlex, ground=QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner, Nemo, QQ, katsura 8",
        result=compute_gb(Groebner.katsuran(8, ordering=:degrevlex, ground=Nemo.QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, QQ, eco 10",
        result=compute_gb(Groebner.eco10(ordering=:degrevlex, ground=QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, QQ, cyclic 7",
        result=compute_gb(Groebner.cyclicn(7, ordering=:degrevlex, ground=QQ))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, QQ, noon 8",
        result=compute_gb(Groebner.noonn(8, ordering=:degrevlex, ground=QQ), 3)
    )
)

#! format: off
R,(t1,t2,t3,a,b,c) = polynomial_ring(QQ, ["t1","t2","t3","a", "b", "c"], ordering=:degrevlex)
hexapod = [1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1065102000*a^2*t1-1566200000*a^2*t2+359610000*a^2*t3-4000000*a*b*t2-1574352000*a*b*t3+4000000*a*c*t1+273640000*a*c*t3-1065102000*b^2*t1+8152000*b^2*t2+355610000*b^2*t3-1574352000*b*c*t1-273640000*b*c*t2-791462000*c^2*t1-1566200000*c^2*t2+355610000*c^2*t3+740236705137*a^2-279943961360*a*b+47071636200*a*c+1574352000*a*t1-273640000*a*t2+126292488913*b^2+837307375312*b*c+4000000*b*t1-273640000*b*t3+612513941897*c^2+4000000*c*t2-1574352000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-624135247952*a-50784764200*b-283060057360*c-791462000*t1+8152000*t2+359610000*t3+165673, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2-1889130000*a^2*t1-139016000*a^2*t2+357608000*a^2*t3+550492000*a*b*t3+1500376000*a*c*t3-1889130000*b^2*t1-689508000*b^2*t2+357608000*b^2*t3+550492000*b*c*t1-1500376000*b*c*t2-388754000*c^2*t1-139016000*c^2*t2+357608000*c^2*t3+740396599024*a^2+98430171568*a*b+268273230304*a*c-550492000*a*t1-1500376000*a*t2+854420557476*b^2-2714848476*b*c-1500376000*b*t3-114024022072*c^2+550492000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2+624263610988*a-268273230304*b+98430171568*c-388754000*t1-689508000*t2+357608000*t3-63620, 4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2-3295636000*a^2*t1+6825304000*a^2*t2+1438448000*a^2*t3-16000000*a*b*t2+4096192000*a*b*t3+16000000*a*c*t1+4906624000*a*c*t3-3295636000*b^2*t1+2729112000*b^2*t2+1422448000*b^2*t3+4096192000*b*c*t1-4906624000*b*c*t2+1610988000*c^2*t1+6825304000*c^2*t2+1422448000*c^2*t3+2962666483625*a^2+722869290752*a*b+875649162944*a*c-4096192000*a*t1-4906624000*a*t2+513760438633*b^2-3361285532000*b*c+16000000*b*t1-4906624000*b*t3+2443184693353*c^2+16000000*c*t2+4096192000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2-2498705324448*a-879018458944*b+741978122752*c+1610988000*t1+2729112000*t2+1438448000*t3+440361,4000000*a^2*t1^2+4000000*a^2*t2^2+4000000*a^2*t3^2+4000000*b^2*t1^2+4000000*b^2*t2^2+4000000*b^2*t3^2+4000000*c^2*t1^2+4000000*c^2*t2^2+4000000*c^2*t3^2+3295636000*a^2*t1+6824896000*a^2*t2+1430432000*a^2*t3+4094592000*a*b*t3-4906624000*a*c*t3+3295636000*b^2*t1+2730304000*b^2*t2+1430432000*b^2*t3+4094592000*b*c*t1+4906624000*b*c*t2-1610988000*c^2*t1+6824896000*c^2*t2+1430432000*c^2*t3+2961910911797*a^2+732129427968*a*b-877323997696*a*c-4094592000*a*t1+4906624000*a*t2+516620569397*b^2+3361357491776*b*c+4906624000*b*t3+2445290017525*c^2+4094592000*c*t3+4000000*t1^2+4000000*t2^2+4000000*t3^2+2499114213824*a+877323997696*b+732129427968*c-1610988000*t1+2730304000*t2+1430432000*t3-324875, 1000000*a^2*t1^2+1000000*a^2*t2^2+1000000*a^2*t3^2+1000000*b^2*t1^2+1000000*b^2*t2^2+1000000*b^2*t3^2+1000000*c^2*t1^2+1000000*c^2*t2^2+1000000*c^2*t3^2+1889602000*a^2*t1-138926000*a^2*t2+359604000*a^2*t3-4000000*a*b*t2+550036000*a*b*t3+4000000*a*c*t1-1500228000*a*c*t3+1889602000*b^2*t1-688962000*b^2*t2+355604000*b^2*t3+550036000*b*c*t1+1500228000*b*c*t2+389374000*c^2*t1-138926000*c^2*t2+355604000*c^2*t3+740903906549*a^2+99175424872*a*b-265964790856*a*c-550036000*a*t1+1500228000*a*t2+854030749541*b^2+2874521168*b*c+4000000*b*t1+1500228000*b*t3-114557203083*c^2+4000000*c*t2+550036000*c*t3+1000000*t1^2+1000000*t2^2+1000000*t3^2-623884900400*a+270522742856*b+97519648872*c+389374000*t1-688962000*t2+359604000*t3+55909, 250000*a^2*t1^2+250000*a^2*t2^2+250000*a^2*t3^2+250000*b^2*t1^2+250000*b^2*t2^2+250000*b^2*t3^2+250000*c^2*t1^2+250000*c^2*t2^2+250000*c^2*t3^2+266341000*a^2*t1-391502000*a^2*t2+89402000*a^2*t3-393620000*a*b*t3-68228000*a*c*t3+266341000*b^2*t1+2118000*b^2*t2+89402000*b^2*t3-393620000*b*c*t1+68228000*b*c*t2+198113000*c^2*t1-391502000*c^2*t2+89402000*c^2*t3+184958257568*a^2-70380830480*a*b-12199439312*a*c+393620000*a*t1+68228000*a*t2+31688927488*b^2-209385275032*b*c+68228000*b*t3+153269490056*c^2-393620000*c*t3+250000*t1^2+250000*t2^2+250000*t3^2+156251491928*a+12199439312*b-70380830480*c+198113000*t1+2118000*t2+89402000*t3+159976]
#! format: on
push!(suite, (problem_name="groebner, AA, QQ, hexapod", result=compute_gb(hexapod, 2)))

function multimodular_gb_problem(nbits; np=AbstractAlgebra)
    R, (x1, x2, x3, x4) =
        polynomial_ring(np.QQ, ["x1", "x2", "x3", "x4"], ordering=:degrevlex)
    nbits_per_prime = 31
    nprimes = max(div(nbits, nbits_per_prime), 1)
    N = prod(map(BigInt, Primes.nextprimes(2^31 - 100, nprimes)))
    @info "Constructing a multi-modular problem" np nbits nbits_per_prime nprimes
    system = [
        x1 + x2 + x3 + x4,
        x1 * x2 + x1 * x3 + x1 * x4 + x2 * x3 + x2 * x4 + x3 * x4,
        x1 * x2 * x3 + x1 * x2 * x4 + x1 * x3 * x4 + x2 * x3 * x4,
        x1 * x2 * x3 * x4 + N
    ]
    system
end

push!(
    suite,
    (
        problem_name="groebner, AA, QQ, 1000-bit output coeffs",
        result=compute_gb(multimodular_gb_problem(1000, np=AbstractAlgebra))
    )
)
push!(
    suite,
    (
        problem_name="groebner, AA, QQ, 10000-bit output coeffs",
        result=compute_gb(multimodular_gb_problem(10000, np=AbstractAlgebra))
    )
)
push!(
    suite,
    (
        problem_name="groebner, Nemo, QQ, 1000-bit output coeffs",
        result=compute_gb(multimodular_gb_problem(1000, np=Nemo))
    )
)
push!(
    suite,
    (
        problem_name="groebner, Nemo, QQ, 10000-bit output coeffs",
        result=compute_gb(multimodular_gb_problem(10000, np=Nemo))
    )
)

function compute_normalforms(system)
    R = parent(system[1])
    gb = Groebner.groebner(system)
    times = []
    trials = 7
    for _ in 1:trials
        GC.gc()
        time = @elapsed begin
            n1 = normalform(gb, system)
            n2 = normalform(gb, gb)
        end
        push!(times, time)
    end
    minimum(times)
end

# Compute normal forms over integers modulo a prime
push!(
    suite,
    (
        problem_name="normalform, AA, GF(2^31-1), cyclic 7",
        result=compute_normalforms(
            Groebner.cyclicn(7, ordering=:degrevlex, ground=GF(2^31 - 1))
        )
    )
)
push!(
    suite,
    (
        problem_name="normalform, AA, GF(103), cyclic 8",
        result=compute_normalforms(
            Groebner.cyclicn(8, ordering=:degrevlex, ground=GF(103))
        )
    )
)
push!(
    suite,
    (
        problem_name="normalform, Nemo, GF(103), cyclic 8",
        result=compute_normalforms(
            Groebner.cyclicn(8, ordering=:degrevlex, ground=Nemo.GF(103))
        )
    )
)
push!(
    suite,
    (
        problem_name="normalform, AA, QQ, katsura 9",
        result=compute_normalforms(Groebner.katsuran(9, ordering=:degrevlex, ground=QQ))
    )
)

function dump_results(file, key)
    open(file, "w") do out
        println(out, key)
        for problem in suite
            problem_name = problem.problem_name
            result = problem.result
            println(out, "$problem_name, $result")
        end
    end
end
