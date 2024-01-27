## Benchmark results

2024-01-27T02:46:22.712

Benchmarked backend: maple

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 60 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.97|
|chandra 12|4.18|
|chandra 13|25.59|
|chandra 14| - |
|cyclic 7|0.30|
|cyclic 8|1.66|
|cyclic 9| - |
|dummy|0.18|
|eco 11|0.54|
|eco 12|2.35|
|eco 13|20.51|
|eco 14| - |
|henrion 5|0.19|
|henrion 6|0.25|
|henrion 7|7.20|
|henrion 8| - |
|katsura 10|1.73|
|katsura 11|15.16|
|katsura 12| - |
|katsura 13| - |
|noon 7| - |
|noon 8|2.09|
|noon 9|19.33|
|noon 10| - |
|reimer 6|0.23|
|reimer 7|1.75|
|reimer 8|54.01|

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
