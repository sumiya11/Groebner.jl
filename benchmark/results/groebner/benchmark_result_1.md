## Benchmark results

2024-01-27T02:44:16.764

Benchmarked backend: groebner

Benchmark suite: Integers modulo 2^30 + 3

- Workers: 16
- Timeout: 60 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 11|0.76|
|chandra 12|3.69|
|chandra 13| - |
|chandra 14| - |
|cyclic 7|0.11|
|cyclic 8|1.62|
|cyclic 9| - |
|dummy|0.00|
|eco 11|0.33|
|eco 12|2.61|
|eco 13|17.10|
|eco 14| - |
|henrion 5|0.00|
|henrion 6|0.03|
|henrion 7|2.70|
|henrion 8| - |
|katsura 10|0.86|
|katsura 11| - |
|katsura 12| - |
|katsura 13| - |
|noon 7|0.18|
|noon 8|1.55|
|noon 9| - |
|noon 10| - |
|reimer 6|0.07|
|reimer 7|0.80|
|reimer 8| - |

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
