## Benchmark results

2024-01-27T07:56:37.996

Benchmarked backend: maple

Benchmark suite: The rationals

- Workers: 8
- Timeout: 3600 s
- Aggregated over: 3 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|5.53|
|chandra 10|25.24|
|chandra 11|137.17|
|chandra 12|927.56|
|chandra 13| - |
|cyclic 7|1.60|
|cyclic 8|22.57|
|cyclic 9| - |
|dummy|0.16|
|eco 10|0.78|
|eco 11|5.30|
|eco 12|39.84|
|eco 13|420.79|
|henrion 6|1.91|
|henrion 7|3205.31|
|ipp|82.50|
|katsura 9|8.71|
|katsura 10|93.27|
|katsura 11|1439.72|
|noon 7|1.18|
|noon 8|9.32|
|noon 9|75.82|
|reimer 6|0.80|
|reimer 7|31.54|
|reimer 8|3425.27|

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
