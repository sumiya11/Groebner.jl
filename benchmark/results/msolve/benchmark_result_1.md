## Benchmark results

2024-01-27T02:48:22.410

Benchmarked backend: msolve

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 60 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.96|
|chandra 12|4.78|
|chandra 13|30.71|
|chandra 14| - |
|cyclic 7|0.15|
|cyclic 8|1.48|
|cyclic 9| - |
|dummy|0.01|
|eco 11|0.58|
|eco 12|3.92|
|eco 13|23.67|
|eco 14| - |
|henrion 5|0.02|
|henrion 6|0.08|
|henrion 7|3.91|
|henrion 8| - |
|katsura 10|1.19|
|katsura 11|10.34|
|katsura 12| - |
|katsura 13| - |
|noon 7|0.23|
|noon 8|2.05|
|noon 9|25.39|
|noon 10| - |
|reimer 6| - |
|reimer 7|1.32|
|reimer 8|47.06|

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
