
import AbstractAlgebra

using Base.Threads

import Pkg
using CpuId
using Printf
using BenchmarkTools
include("../../src/Groebner.jl")
include("../systems/for_gleb/parser.jl")

using Logging
global_logger(ConsoleLogger(stderr, Logging.Error))

BENCHMARK_SAMPLES = 4
GROEBNER_PARAMS = (linalg=:prob, )

ground = AbstractAlgebra.QQ

systems = [
    ("siwr", read_SIWR()),
    ("seaijrc", read_SEAIJRC()),
    ("cyclic 7",  Groebner.cyclicn(7, ground=ground)), 
    ("cyclic 8",  Groebner.cyclicn(8, ground=ground)), 
    # ("cyclic 9",  Groebner.cyclicn(9, ground=ground)), 
    # ("katsura 10",Groebner.katsuran(10, ground=ground)), 
    ("katsura 9",Groebner.katsuran(9, ground=ground)), 
    ("katsura 10",Groebner.katsuran(10, ground=ground)), 
    ("katsura 11",Groebner.katsuran(11, ground=ground)), 
    ("eco 11",    Groebner.eco11(ground=ground)),         
    ("eco 12",    Groebner.eco12(ground=ground)),         
    ("eco 13",    Groebner.eco13(ground=ground)),         
    # ("noon 7",    Groebner.noonn(7, ground=ground)),  
    ("noon 8",    Groebner.noonn(8, ground=ground)),  
    ("noon 9",    Groebner.noonn(9, ground=ground)),   
    # ("henrion 5", Groebner.henrion5(ground=ground)),  
    ("henrion 6", Groebner.henrion6(ground=ground)),  
    # ("henrion 7", Groebner.henrion7(ground=ground)),  
    ("reimer 6",  Groebner.reimern(6, ground=ground)), 
    ("reimer 7",  Groebner.reimern(7, ground=ground)), 
    # ("reimer 8",  Groebner.reimern(8, ground=ground)), 
]

function benchmark_system_groebner(system)
    system = Groebner.change_ordering(system, :degrevlex)
    # Reduce compilation times
    global GROEBNER_PARAMS
    global BENCHMARK_SAMPLES
    println("Basis size = ", length(Groebner.groebner(system; GROEBNER_PARAMS...)))
    # run benchmarks
    bench = @benchmarkable Groebner.groebner($system; $GROEBNER_PARAMS...) samples=BENCHMARK_SAMPLES
    median(run(bench))
end

function run_system_groebner(name, system, resulting_md)
    x = benchmark_system_groebner(system)
    println("---------------\n$name")
    println("groebner")
    println(x)
    Groebner.printall()
    println("---------------")
    resulting_md *= string(x) * "|"
end

#= 
    A part of this benchmarking script is taken 
    from StructuralIdentifiability.jl benchmarks
    https://github.com/SciML/StructuralIdentifiability.jl/blob/master/benchmarking/run_benchmarks.jl
=#

function run_f4_ff_degrevlex_benchmarks()

    resulting_md = ""
    resulting_md *= "|System|Groebner.jl|\n"

    for (name, system) in systems
        resulting_md *= "|$name|"
        run_system_groebner(name, system, resulting_md)
        resulting_md *= "\n"
    end

    global GROEBNER_PARAMS
    resulting_md *= "\n*Groebner parameters:*\n"
    resulting_md *= string(GROEBNER_PARAMS)
    resulting_md *= "\n*Benchmarking environment:*\n\n"
    resulting_md *= "* Total RAM (Mb): $(Sys.total_memory() / 2^20)\n"
    resulting_md *= "* Processor: $(cpubrand())\n"
    resulting_md *= "* Julia version: $(VERSION)\n\n"
    resulting_md *= "Versions of the dependencies:\n\n"

    deps = Pkg.dependencies()
    stid_info = deps[findfirst(x -> x.name == "Groebner", deps)]
    for (s, uid) in stid_info.dependencies
        if deps[uid].version !== nothing
            resulting_md *= "* $s : $(deps[uid].version)\n"
        end
    end

    open("groebner_benchmark_result.md", "w") do io
        write(io, resulting_md)
    end
end

run_f4_ff_degrevlex_benchmarks()

#=

Basis size = 8
---------------
siwr
groebner
TrialEstimate(2.634 s)
run success,
PRIMES = 36
F4TIME = 0.810926650786335
RECTIME = 0.0010346620828865004
CORRTIME = 0.1880386871307786
---------------
Basis size = 10
---------------
seaijrc
groebner
TrialEstimate(8.680 s)
run success,
PRIMES = 36
F4TIME = 0.6017336825580868
RECTIME = 0.0007489715082925581
CORRTIME = 0.3975173459336206
---------------
Basis size = 209
---------------
cyclic 7
groebner
TrialEstimate(4.930 s)
run success,
PRIMES = 36
F4TIME = 0.7025681972591288
RECTIME = 0.28426320117101234
CORRTIME = 0.013168601569858909
---------------
Basis size = 372
---------------
cyclic 8
groebner
TrialEstimate(142.003 s)
run success,
PRIMES = 68
F4TIME = 0.9114366086922844
RECTIME = 0.08655601087181362
CORRTIME = 0.0020073804359020068
---------------
Basis size = 272
---------------
katsura 9
groebner
TrialEstimate(10.409 s)
run success,
PRIMES = 36
F4TIME = 0.5526910379889678
RECTIME = 0.4253887070744696
CORRTIME = 0.02192025493656265
---------------
Basis size = 537
---------------
katsura 10
groebner
TrialEstimate(106.181 s)
run success,
PRIMES = 68
F4TIME = 0.5905210911487123
RECTIME = 0.3998071079679621
CORRTIME = 0.009671800883325641
---------------
Basis size = 1050
---------------
katsura 11
groebner
TrialEstimate(1446.875 s)
run success,
PRIMES = 132
F4TIME = 0.6123534590787052
RECTIME = 0.384348475953135
CORRTIME = 0.0032980649681598495
---------------
Basis size = 389
---------------
eco 11
groebner
TrialEstimate(9.192 s)
run success,
PRIMES = 20
F4TIME = 0.6732230923728532
RECTIME = 0.3011675205799542
CORRTIME = 0.025609387047192483
---------------
Basis size = 743
---------------
eco 12
groebner
TrialEstimate(47.515 s)
run success,
PRIMES = 20
F4TIME = 0.8316717141074981
RECTIME = 0.15386372995807887
CORRTIME = 0.014464555934423026
---------------
Basis size = 1465
---------------
eco 13
groebner
TrialEstimate(393.127 s)
run success,
PRIMES = 36
F4TIME = 0.8564150884215842
RECTIME = 0.1370970912105766
CORRTIME = 0.006487820367839246
---------------
Basis size = 1338
---------------
noon 8
groebner
TrialEstimate(6.172 s)
run success,
PRIMES = 4
F4TIME = 0.6901164045660294
RECTIME = 0.24293745778454434
CORRTIME = 0.06694613764942624
---------------
Basis size = 3682
---------------
noon 9
groebner
TrialEstimate(61.782 s)
run success,
PRIMES = 5
F4TIME = 0.8182580732755192
RECTIME = 0.15482847986346904
CORRTIME = 0.026913446861011624
---------------
Basis size = 90
---------------
henrion 6
groebner
TrialEstimate(2.413 s)
run success,
PRIMES = 36
F4TIME = 0.4710110842875254
RECTIME = 0.5050115446227048
CORRTIME = 0.023977371089769797
---------------
Basis size = 95
---------------
reimer 6
groebner
TrialEstimate(848.489 ms)
run success,
PRIMES = 12
F4TIME = 0.5978586617686001
RECTIME = 0.3580209335948437
CORRTIME = 0.04412040463655626
---------------
Basis size = 227
---------------
reimer 7
groebner
TrialEstimate(34.105 s)
run success,
PRIMES = 36
F4TIME = 0.7555069824554617
RECTIME = 0.23191510291276185
CORRTIME = 0.012577914631776422
---------------

=#
