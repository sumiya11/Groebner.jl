## Benchmark results

2024-01-12T12:45:58.856

Benchmarked backend: maple

Benchmark suite: The rationals

- Workers: 8
- Timeout: 2000 s
- Aggregated over: 1 runs

**All timings in seconds.**

|Model|Total, s|
|:----|---|
|chandra 9|5.85|
|chandra 10|27.81|
|chandra 11|151.60|
|chandra 12|1166.48|
|chandra 13| - |
|cyclic 7|1.50|
|cyclic 8|24.92|
|cyclic 9| - |
|dummy|0.01|
|eco 10|0.72|
|eco 11|4.72|
|eco 12|39.68|
|eco 13|532.73|
|henrion 6|1.96|
|henrion 7| - |
|ipp|196.19|
|katsura 9|9.07|
|katsura 10|102.40|
|katsura 11|1479.51|
|noon 7|1.02|
|noon 8|9.62|
|noon 9|102.74|
|reimer 6|0.65|
|reimer 7|37.55|
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
* AbstractAlgebra : 0.34.7
* Nemo : 0.38.3
* Atomix : 0.1.0
* ExprTools : 0.1.10
* PrettyTables : 2.3.1
