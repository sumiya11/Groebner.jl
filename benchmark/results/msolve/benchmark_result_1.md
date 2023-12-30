## Benchmark results

2023-12-30T23:39:47.549

Benchmarked backend: msolve

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 20 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 7| - |
|cyclic 8| - |
|cyclic 9| - |
|cyclic 10| - |
|dummy| - |
|eco 11| - |
|eco 12| - |
|eco 13| - |
|eco 14| - |
|henrion 5| - |
|henrion 6| - |
|henrion 7| - |
|katsura 10| - |
|katsura 11| - |
|katsura 12| - |
|katsura 13| - |
|noon 7| - |
|noon 8| - |
|noon 9| - |
|noon 10| - |
|reimer 6| - |
|reimer 7| - |
|reimer 8| - |
|reimer 9| - |

*Benchmarking environment:*

* Total RAM (GiB): 188
* Processor: Intel(R) Xeon(R) Gold 6130 CPU @ 2.10GHz
* Julia version: 1.9.2

Versions of the dependencies:

* Primes : 0.5.5
* TimerOutputs : 0.5.23
* PrecompileTools : 1.2.0
* MultivariatePolynomials : 0.5.3
* Combinatorics : 1.0.2
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
