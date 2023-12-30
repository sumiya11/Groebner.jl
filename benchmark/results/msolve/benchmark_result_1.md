## Benchmark results

2023-12-30T12:48:59.245

Benchmarked backend: msolve
Benchmark suite: Integers modulo 2^30 + 3

- Workers: 4
- Timeout: 1200 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total|
|-----|---|
|cyclic 10| - |
|cyclic 7|0.14|
|cyclic 8|1.49|
|cyclic 9|268.64|
|dummy|0.02|
|eco 11|0.53|
|eco 12|2.96|
|eco 13|24.89|
|eco 14|273.12|
|henrion 5|0.02|
|henrion 6|0.08|
|henrion 7|4.02|
|katsura 10|1.11|
|katsura 11|11.05|
|katsura 12|110.44|
|katsura 13|823.39|
|noon 10|318.72|
|noon 7|0.23|
|noon 8|2.18|
|noon 9|24.19|
|reimer 6|0.08|
|reimer 7|1.23|
|reimer 8|48.93|
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