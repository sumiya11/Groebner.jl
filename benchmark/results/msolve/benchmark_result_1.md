## Benchmark results

2024-02-03T15:43:02.646

Benchmarked backend: msolve

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.90|
|chandra 12|4.39|
|chandra 13|22.34|
|chandra 14|185.78|
|cyclic 7|0.15|
|cyclic 8|1.38|
|cyclic 9|125.09|
|cyclic 10| - |
|dummy|0.01|
|eco 11|0.54|
|eco 12|2.98|
|eco 13|18.39|
|eco 14|147.85|
|henrion 5|0.02|
|henrion 6|0.08|
|henrion 7|3.52|
|henrion 8|907.76|
|katsura 10|1.18|
|katsura 11|7.47|
|katsura 12|70.53|
|katsura 13|483.17|
|noon 7|0.23|
|noon 8|2.35|
|noon 9|20.09|
|noon 10|226.64|
|noon 11| - |
|reimer 6|0.07|
|reimer 7|1.30|
|reimer 8|42.89|
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
* HostCPUFeatures : 0.1.16
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
